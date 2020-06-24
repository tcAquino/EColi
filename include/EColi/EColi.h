//
// EColi.h
// EColi
//
// Created by Tomas Aquino on 8/2/19.
// Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Classes to build CTRW model for EColi dynamics
// Region values correspond to
//    0 : water column
//    1 : hyporheic zone
//    2 : bed
//    3 : Absorbed (past maximum distance)

#ifndef EColi_h
#define EColi_h

#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <utility>
#include <stdexcept>
#include "general/useful.h"
#include "Stochastic/Random.h"

namespace ecoli
{
  // EColi state for CTRW
  struct State_EColi
  {
    State_EColi
    (double position, double mass,
     std::size_t region, double time)
    : position{ position }
    , mass{ mass }
    , region{ region }
    , time{ time }
    {}

    double position;     // downstream position
    double mass;         // amount of mass carried by particle
    std::size_t region;  // 0: Water column, 1: hyporheic zone, 2: bed, 3: absorbed
    double time;         // Each state has its own time, corresponding to its onset
  };

  // First order decay with region-dependent rates
  class Reaction_EColi_Decay
  {
  public:
    const struct Parameters
    {
      Parameters(double rate_w, double rate_h, double rate_b)
      : rate( 3 )
      {
        rate[0] = rate_w;
        rate[1] = rate_h;
        rate[2] = rate_b;
      }

      std::vector< double > rate;
    }
    parameters;

    Reaction_EColi_Decay(Parameters parameters)
    : parameters{ parameters }
    {}

    template <typename State_t>
    void operator() (State_t& state, double exposure_time) const
    {
      state.mass *= std::exp(-parameters.rate[state.region] * exposure_time);
    }
  };

  // Advection class to handle 1d piecewise advection
  // E.g. a velocity time series v(t) = { {t1, v1}, {t2,v2} }
  // corresponds to v(t)=v1 for t1<=t<t2, v(t)=v2 for t>t2.
  class Advection_EColi_Piecewise
  {
  public:
    using Container = std::vector<std::tuple<double, double, bool>>;
    using Element = Container::value_type;
    using Iterator = Container::const_iterator;

    Advection_EColi_Piecewise(Container advection)
    : advection{ advection }
    {}

    Advection_EColi_Piecewise(std::string const& filename)
    : advection{ read(filename) }
    {}
    
    Advection_EColi_Piecewise
    (std::string const& filename, std::string const& filename_surge)
    : advection{ read(filename, filename_surge) }
    {}

    auto cbegin() const
    { return advection.cbegin(); }

    auto cend() const
    { return advection.cend(); }

    // iterator to (onset time,velocity) corresponding to given time
    auto iterator(double time) const
    {
      return --std::upper_bound
      (advection.cbegin(), advection.cend(), time,
       [](double const& val, Element const& a2)
       { return val < std::get<0>(a2); });
    }

    // (onset time,velocity,surge) corresponding to given time
    Element operator()(double time) const
    { return *iterator(time); }

    // onset time corresponding to given time
    double time(double time) const
    { return std::get<0>(*iterator(time)); }

    // velocity corresponding to given time
    double velocity(double time) const
    { return std::get<1>(*iterator(time)); }
    
    // surge corresponding to given time
    double surge(double time) const
    { return std::get<2>(*iterator(time)); }
    
    // (onset time,velocity,surge) corresponding to given iterator
    Element operator()(Iterator it) const
    { return *it; }

    // surge corresponding to given iterator
    double time(Iterator it) const
    { return std::get<0>(*it); }

    // surge corresponding to given iterator
    double velocity(Iterator it) const
    { return std::get<1>(*it); }
    
    // surge corresponding to given time
    double surge(Iterator it) const
    { return std::get<2>(*it); }

    // Jump length corresponding to moving according to v(t)
    // for t in [state.time, state.time+exposure_time],
    // \int_state.time^{state.time+exposure_time} dt v(t)
    template <typename State_t>
    double advection_jump
    (State_t const& state, double exposure_time) const
    {
      // Advection happens only in the water column
      if (state.region == 0)
      {
        auto const it_min = iterator(state.time);
        auto const it_max = iterator(state.time + exposure_time);

        // Handle final time being within the first velocity time window
        if (it_min == it_max)
          return velocity(it_min)*exposure_time;

        // Otherwise, compute jump to end of first velocity time window
        double delta_t = time(it_min+1)-state.time;
        double jump = velocity(it_min)*delta_t;

        // Intermediate velocity time windows
        for (auto it = it_min + 1; it < it_max; ++it)
        {
          delta_t = time(it+1)-time(it);
          jump += velocity(it)*delta_t;
        }

        // Final velocity window
        delta_t = state.time+exposure_time-time(it_max);
        jump += velocity(it_max)*delta_t;

        return jump;
      }
      else
      {
        return 0.;
      }
    }

  private:
    const Container advection;

    Container read(std::string const& filename)
    {
      std::ifstream file{ filename };
      if (!file.is_open())
        throw useful::open_read_error(filename);
      Container advection;
      double time, val;
      while(file >> time && file >> val)
        advection.push_back({ time, val, 0 });
      file.close();

      return advection;
    }
    
    Container read
    (std::string const& filename, std::string const& filename_surge)
    {
      std::ifstream file{ filename };
      if (!file.is_open())
        throw useful::open_read_error(filename);
      std::ifstream file_surge{ filename_surge };
      if (!file_surge.is_open())
        throw useful::open_read_error(filename_surge);
      Container advection;
      double time, val;
      bool surge;
      while(file >> time && file >> val && file_surge >> surge)
        advection.push_back({ time, val, surge });
      file.close();
      file_surge.close();

      return advection;
    }
  };

  // Transitions helper class to be used with Transitions class below
  // Computes residence time and region to transition to
  // Handles piecewise advection in water column,
  // exponential residence time in hyporheic zone due to dispersion/advection,
  // and stable residence time bed
  template
  <typename Time_WaterColumn, typename Time_Hyporheic>
  class TransitionHelper_EColi_PiecewiseAdvection
  {
  public:
    // Constructors
    const struct Parameters
    {
      double time_h_disp;                 // Mean dispersive residence time in hyporheic zone
      double time_b;                      // Characteristic residence time in bed
      double alpha_b;                     // Exponent of residence time in bed
      Time_WaterColumn time_watercolumn;  // Full mean residence time in water column
      Time_Hyporheic time_hyporheic_adv;  // Mean residence time by advection in hyporheic zone
      double prob_resusp{ 0.5 };          // Probability of transition to water column (dispersion)
      double prob_resusp_adv{ 1. };       // Probability of transition to water column (advection)
    } parameters;

    TransitionHelper_EColi_PiecewiseAdvection(Parameters parameters)
    : parameters{ parameters }
    {}

    // Return the time to the next transition and the new region as a pair
    template <typename State_t, typename Advection>
    auto operator()
    (State_t const& state, Advection const& advection)
    {
      switch (state.region)
      {
        case 0:
          return transition_watercolumn(state,advection);

        case 1:
          return transition_hyporheic(state,advection);

        case 2:
          return transition_bed(state,advection);

        default:
          throw std::runtime_error{ "Invalid region." };
      }
    }

  private:
    // Random number generation
    std::mt19937 rng{ std::random_device{}() };
    std::exponential_distribution<double> dist_exp{};
    std::bernoulli_distribution dist_ber_disp{
      1.-parameters.prob_resusp };
    std::bernoulli_distribution dist_ber_adv{
      1.-parameters.prob_resusp_adv };
    stochastic::skewedlevystable_distribution<double>
      dist_lev{ parameters.alpha_b };

    double time_watercolumn(double advection, bool surge)
    { return parameters.time_watercolumn(advection, surge)*
        dist_exp(rng); }

    double time_hyporheic_adv(double advection, bool surge)
    { return parameters.time_hyporheic_adv(advection, surge)*
        dist_exp(rng); }

    double time_bed()
    { return parameters.time_b*dist_lev(rng); }

    double time_hyporheic_disp()
    { return parameters.time_h_disp*dist_exp(rng); }

    std::size_t region_hyporheic_disp()
    { return 2*dist_ber_disp(rng); }

    std::size_t region_hyporheic_adv(bool surge)
    { return surge ? 0 : 2*dist_ber_adv(rng); }

    template <typename State_t, typename Advection>
    auto transition_watercolumn
    (State_t const& state, Advection const& advection)
    {
      std::pair<std::size_t, double> region_time;

      // Transition to hyporheic zone
      region_time.first = 1;

      // Advection window for current time and last window
      auto const adv_it_min = advection.iterator(state.time);
      auto const adv_it_max = advection.cend()-1;

      // If first window is the same as last window
      // return first advection time
      double time_advection_window = time_watercolumn(
        advection.velocity(adv_it_min), advection.surge(adv_it_min));
      if (adv_it_min == adv_it_max)
      {
        region_time.second = time_advection_window;
        return region_time;
      }

      // Check for resuspension by advection until end of first window
      double delta_t = advection.time(adv_it_min+1)-state.time;
      if (time_advection_window < delta_t)
      {
        region_time.second = time_advection_window;
        return region_time;
      }
      double time_advection = delta_t;

      // Check for transition in intermediate velocity windows
      for (auto it = adv_it_min + 1; it < adv_it_max; ++it)
      {
        time_advection_window = time_watercolumn(
          advection.velocity(it), advection.surge(it));
        delta_t = advection.time(it+1) - advection.time(it);
        if (time_advection_window < delta_t)
        {
          region_time.second = time_advection + time_advection_window;
          return region_time;
        }
        time_advection += delta_t;
      }

      // Same for last velocity window
      region_time.second = time_advection+
        time_watercolumn(advection.velocity(adv_it_max),
          advection.surge(adv_it_max));
      return region_time;
    }

    // If time by advection < than by dispersion, transition probabilities to water column/bed by advection
    // Otherwise, transition probabilities to water column/bed by dispersion
    template <typename State_t, typename Advection>
    auto transition_hyporheic
    (State_t const& state, Advection const& advection)
    {
      std::pair<std::size_t, double> region_time;

      // Time to transition by dispersion
      const double time_disp = time_hyporheic_disp();

      // Advection window for current time and for the time corresponding to the dispersion transition
      auto const adv_it_min = advection.iterator(state.time);
      auto const adv_it_max = advection.iterator(state.time+time_disp);

      // Compute time to resuspension by advection for first advection window
      double time_adv_window = time_hyporheic_adv(
        advection.velocity(adv_it_min), advection.surge(adv_it_min));
      bool surge = advection.surge(adv_it_min);

      // If first window is the same as last window,
      // find next region and time and return
      if (adv_it_min == adv_it_max)
      {
        //	 Transition by advection
        if (time_adv_window <= time_disp)
        {
          region_time.first = region_hyporheic_adv(surge);
          region_time.second = time_adv_window;
        }
        //	 Transition by dispersion
        else
        {
          region_time.first = region_hyporheic_disp();
          region_time.second = time_disp;
        }
        return region_time;
      }

      // Check for resuspension by advection until end of first window
      double delta_t = advection.time(adv_it_min+1)-state.time;
      if (time_adv_window <= delta_t)
      {
        region_time.first = region_hyporheic_adv(
          advection.surge(adv_it_min));
        region_time.second = time_adv_window;
        return region_time;
      }
      double time_advection = delta_t;

      // Same for intermediate velocity windows
      for (auto it = adv_it_min + 1; it < adv_it_max; ++it)
      {
        surge = advection.surge(it);
        time_adv_window = time_hyporheic_adv(
          advection.velocity(it), surge);
        delta_t = advection.time(it+1)-advection.time(it);
        if (time_adv_window <= delta_t)
        {
          region_time.first = region_hyporheic_adv(surge);
          region_time.second = time_advection + time_adv_window;
          return region_time;
        }
        time_advection += delta_t;
      }

      // Same for last velocity window
      surge = advection.surge(adv_it_max);
      time_adv_window = time_hyporheic_adv(
        advection.velocity(adv_it_max), surge);
      if (time_advection + time_adv_window <= time_disp)
      {
        region_time.first = region_hyporheic_adv(surge);
        region_time.second = time_advection + time_adv_window;
      }
      else
      {
        region_time.first = region_hyporheic_disp();
        region_time.second = time_disp;
      }
      return region_time;
    }

    template <typename State_t, typename Advection>
    auto transition_bed
    (State_t const& state, Advection const& advection)
    {
      std::pair<std::size_t, double> region_time;

      // Transition to hyporheic zone
      region_time.first = 1;
      region_time.second = time_bed();

      return region_time;
    }
  };

  // Handle state changes according to the dynamics
  // Piecewise advection,
  // residence times according to TransitionHelper_EColi_PiecewiseAdvection,
  // first order decay with region-dependent rates,
  // and absorption past a given maximum distance
  template <typename TransitionHelper>
  class Transitions_EColi
  {
  private:
    using Reaction = Reaction_EColi_Decay;

  public:
    using Advection = Advection_EColi_Piecewise;
    using Parameters_residence = typename TransitionHelper::Parameters;
    using Parameters_reaction = typename Reaction::Parameters;

    Transitions_EColi
    (Advection advection,
     Parameters_residence residence_time_parameters,
     Parameters_reaction reaction_parameters,
     double max_distance)
    : adv{ advection }
    , transition_helper{ residence_time_parameters }
    , react{ reaction_parameters }
    , max_distance{ max_distance }
    {}

    template <typename State_t>
    void operator() (State_t& state)
    {
      // Ignore absorbed particles
      if (state.region == 3) return;

      // Compute new state properties
      auto region_time = transition_helper(state, adv);
      auto new_region = region_time.first;
      auto exposure_time = region_time.second;
      auto jump = adv.advection_jump(state, exposure_time);

      // Deal with jump going out of bounds
      // or set the new state as computed
      if (!OutOfBounds(state, jump))
      {
        react(state, exposure_time);
        state.region = new_region;
        state.position += jump;
        state.time += exposure_time;
      }
    }

    Advection const& advection()
    { return adv; }

    Reaction const& reaction()
    { return react; }

  private:
    Advection adv;
    TransitionHelper transition_helper;
    Reaction react;
    const double max_distance;

    // Check if jump takes the current position beyond the maximum distance
    // If out of bounds, set the region to 3,
    // set the time to the time of reaching max_distance,
    // and return true.
    // Otherwise, return false
    template <typename State_t>
    bool OutOfBounds(State_t& state, double jump)
    {
      // If out of bounds,
      // find out how long until max_distance was reached, update state, and return true
      if (state.position + jump > max_distance)
      {
        double exposure_time = 0.;
        double pos = state.position;

        // For each advection window not the last recorded
        for (auto adv_it = adv.iterator(state.time);
             adv_it < adv.cend() - 1; ++adv_it)
        {
          // Increase exposure time until the beginning of the next window
          double delta_t = adv.time(adv_it+1)-adv.time(adv_it);
          pos += adv.velocity(adv_it)*delta_t;
          exposure_time += delta_t;

          // If the maximum distance is reached in current window,
          // decrease the exposure time,
          // taking into account the non-traversed portion at the current velocity,
          // set the new state, and return true
          if (pos > max_distance)
          {
            exposure_time -= (pos-max_distance)/adv.velocity(adv_it);
            react(state, exposure_time);
            state.region = 3;
            state.position = max_distance;
            state.time += exposure_time;
            return 1;
          }
        }

        // If the crossing occurs in the last available advection window,
        // increase time until max_distance with the last velocity and return false
        exposure_time += (max_distance-pos)/adv.velocity(adv.cend()-1);
        react(state, exposure_time);
        state.region = 3;
        state.position = max_distance;
        state.time += exposure_time;
        return 1;
      }
      
      // Otherwise, return false
      return 0;
    }
  };
}

#endif /* EColi_h */

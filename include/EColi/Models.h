//
//  Models.h
//  EColi
//
//  Created by Tomás Aquino on 01/11/2020.
//  Copyright © 2020 Tomas Aquino. All rights reserved.
//

#ifndef Models_EColi_h
#define Models_EColi_h

#include <iomanip>
#include <iostream>
#include <string>
#include "DepositionResuspension.h"
#include "general/useful.h"
#include "Stochastic/Random.h"

namespace ecoli
{
  namespace model_1
  {
    const std::string model_name = { "1" };
    
    struct Parameters
    {
      double time_w_diff;
      double time_h_diff;
      double time_b;
      double alpha_b;
      double advection_coeff_w;
      double advection_coeff_h;
      double prob_resusp;
      double prob_resusp_adv;
      double rate_decay_w;
      double rate_decay_h;
      double rate_decay_b;
      double concentration_inj_w;
      double mass_inj_h;
      double mass_inj_b;
      double time_injection_min;
      double time_injection_max;
      double time_step_injection;
      double max_dist;
      double time_min;
      double time_max;
      std::size_t nr_measurements;
      std::string measure_spacing;
      std::size_t nr_particles;
      
      Parameters(std::string const& filename)
      {
        auto input = useful::open_read(filename);
        useful::read(input, time_w_diff);
        useful::read(input, time_h_diff);
        useful::read(input, time_b);
        useful::read(input, alpha_b);
        useful::read(input, advection_coeff_w);
        useful::read(input, advection_coeff_h);
        useful::read(input, prob_resusp);
        useful::read(input, prob_resusp_adv);
        useful::read(input, rate_decay_w);
        useful::read(input, rate_decay_h);
        useful::read(input, rate_decay_b);
        useful::read(input, concentration_inj_w);
        useful::read(input, mass_inj_h);
        useful::read(input, mass_inj_b);
        useful::read(input, time_injection_min);
        useful::read(input, time_injection_max);
        useful::read(input, time_step_injection);
        useful::read(input, max_dist);
        useful::read(input, time_min);
        useful::read(input, time_max);
        useful::read(input, nr_measurements);
        useful::read(input, measure_spacing);
        useful::read(input, nr_particles);
        input.close();
      }
      
      static void parameter_list()
      {
        std::cout << "Model name: " << model_name << "\n"
                  << "Parameters:\n"
                  << "\ttime_w_diff             : Mean dispersive residence time in water column\n"
                  << "\ttime_h_diff             : Mean dispersive residence time in hyporheic zone\n"
                  << "\ttime_b                  : Characteristic residence time in bed\n"
                  << "\talpha_b                 : Residence time exponent in bed\n"
                  << "\tadvection_coeff_w       : Coefficient of advection-induced rate in water column\n"
                  << "\tadvection_coeff_h       : Coefficient of advection-induced rate in hyporheic zone\n"
                  << "\tprob_resusp             : Probability of transition to water column by diffusion during base flow\n"
                  << "\tprob_resusp_adv         : Probability of transition to water column by advection during base flow\n"
                  << "\trate_decay_w            : Decay rate in water column\n"
                  << "\trate_decay_h            : Decay rate in hyporheic zone\n"
                  << "\trate_decay_b            : Decay rate in bed\n"
                  << "\tconcentration_inj_w     : Concentration (mass if pulse) at upstream boundary in water column\n"
                  << "\tmass_inj_h              : Mass in pulse at upstream boundary in hyporheic zone\n"
                  << "\tmass_inj_b              : Mass in pulse at upstream boundary in bed\n"
                  << "\ttime_injection_min      : Start time of injection\n"
                  << "\ttime_injection_max      : End time of injection (equal to start time for pulse)\n"
                  << "\ttime_step_injection     : Discretization time step for injection (ignored for pulse)\n"
                  << "\tmax_dist                : Maximum downstream distance\n"
                  << "\ttime_min                : Minimum measure time\n"
                  << "\ttime_max                : Maximum measure time\n"
                  << "\tnr_measurements         : Nr of measurement times\n"
                  << "\tmeasure_spacing         : lin - Linear time spacing\n"
                  << "\t                          log - Logarithmic time spacing\n"
                  << "\tnr_particles            : Number of particles per injection\n";
      }
    };
    
    using MeanTime_WaterColumn_Base = MeanTime_Quadratic;
    using MeanTime_Hyporheic_Adv_Base = MeanTime_Linear;
    using HyporeicTransitionRegion_Base = HyporheicTransitionRegion_Probability;
    using HyporeicTransitionRegion_Adv_Base = HyporheicTransitionRegion_Probability;
    
    using MeanTime_WaterColumn_Surge = MeanTime_Infinite;
    using MeanTime_Hyporheic_Adv_Surge = MeanTime_Zero;
    using HyporeicTransitionRegion_Surge = HyporeicTransitionRegion_Base;
    using HyporeicTransitionRegion_Adv_Surge = HyporheicTransitionRegion_Resuspension;
    
    using Time_Bed = stochastic::RNG<stochastic::skewedlevystable_distribution<double>>;

    MeanTime_WaterColumn_Base make_MeanTime_WaterColumn_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_w, 1./parameters.time_w_diff };
    }
    MeanTime_Hyporheic_Adv_Base make_MeanTime_Hyporheic_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_h };
    }
    HyporeicTransitionRegion_Base make_HyporeicTransitionRegion_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp };
    }
    HyporeicTransitionRegion_Adv_Base make_HyporeicTransitionRegion_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp_adv };
    }
    
    MeanTime_WaterColumn_Surge make_MeanTime_WaterColumn_Surge
    (Parameters const& parameters)
    {
      return {};
    }
    MeanTime_Hyporheic_Adv_Surge make_MeanTime_Hyporheic_Adv_Surge
    (Parameters const& parameters)
    {
      return {};
    }
    HyporeicTransitionRegion_Surge make_HyporeicTransitionRegion_Surge
    (Parameters const& parameters)
    {
      return make_HyporeicTransitionRegion_Base(parameters);
    }
    HyporeicTransitionRegion_Adv_Surge make_HyporeicTransitionRegion_Adv_Surge
    (Parameters const& parameters)
    {
      return {};
    }

    Time_Bed make_Time_Bed
    (Parameters const& parameters)
    {
      return Time_Bed::param_type{
        parameters.alpha_b, parameters.time_b };
    }
  }
  
  namespace model_2
  {
    const std::string model_name = { "2" };
    
    struct Parameters
    {
      double time_w_diff;
      double time_h_diff;
      double time_b;
      double alpha_b;
      double advection_coeff_w;
      double advection_coeff_h;
      double prob_resusp;
      double prob_resusp_adv;
      double rate_decay_w;
      double rate_decay_h;
      double rate_decay_b;
      double concentration_inj_w;
      double mass_inj_h;
      double mass_inj_b;
      double time_injection_min;
      double time_injection_max;
      double time_step_injection;
      double max_dist;
      double time_min;
      double time_max;
      std::size_t nr_measurements;
      std::string measure_spacing;
      std::size_t nr_particles;
      
      Parameters(std::string const& filename)
      {
        auto input = useful::open_read(filename);
        useful::read(input, time_w_diff);
        useful::read(input, time_h_diff);
        useful::read(input, time_b);
        useful::read(input, alpha_b);
        useful::read(input, advection_coeff_w);
        useful::read(input, advection_coeff_h);
        useful::read(input, prob_resusp);
        useful::read(input, prob_resusp_adv);
        useful::read(input, rate_decay_w);
        useful::read(input, rate_decay_h);
        useful::read(input, rate_decay_b);
        useful::read(input, concentration_inj_w);
        useful::read(input, mass_inj_h);
        useful::read(input, mass_inj_b);
        useful::read(input, time_injection_min);
        useful::read(input, time_injection_max);
        useful::read(input, time_step_injection);
        useful::read(input, max_dist);
        useful::read(input, time_min);
        useful::read(input, time_max);
        useful::read(input, nr_measurements);
        useful::read(input, measure_spacing);
        useful::read(input, nr_particles);
        input.close();
      }
      
      static void parameter_list()
      {
        std::cout << "Model name: " << model_name << "\n"
                  << "Parameters (default in []):\n"
                  << "\ttime_w_diff             : Mean dispersive residence time in water column\n"
                  << "\ttime_h_diff             : Mean dispersive residence time in hyporheic zone\n"
                  << "\ttime_b                  : Characteristic residence time in bed\n"
                  << "\talpha_b                 : Residence time exponent in bed\n"
                  << "\tadvection_coeff_w       : Coefficient of advection-induced rate in water column\n"
                  << "\tadvection_coeff_h       : Coefficient of advection-induced rate in hyporheic zone\n"
                  << "\tprob_resusp             : Probability of transition to water column by diffusion during base flow\n"
                  << "\tprob_resusp_adv         : Probability of transition to water column by advection during base flow\n"
                  << "\trate_decay_w            : Decay rate in water column\n"
                  << "\trate_decay_h            : Decay rate in hyporheic zone\n"
                  << "\trate_decay_b            : Decay rate in bed\n"
                  << "\tconcentration_inj_w     : Concentration (mass if pulse) at upstream boundary in water column\n"
                  << "\tmass_inj_h              : Mass in pulse at upstream boundary in hyporheic zone\n"
                  << "\tmass_inj_b              : Mass in pulse at upstream boundary in bed\n"
                  << "\ttime_injection_min      : Start time of injection\n"
                  << "\ttime_injection_max      : End time of injection (equal to start time for pulse)\n"
                  << "\ttime_step_injection     : Discretization time step for injection (ignored for pulse)\n"
                  << "\tmax_dist                : Maximum downstream distance\n"
                  << "\ttime_min                : Minimum measure time\n"
                  << "\ttime_max                : Maximum measure time\n"
                  << "\tnr_measurements         : Nr of measurement times\n"
                  << "\tmeasure_spacing         : lin - Linear time spacing\n"
                  << "\t                          log - Logarithmic time spacing\n"
                  << "\tnr_particles            : Number of particles per injection\n";
      }
    };
    
    using MeanTime_WaterColumn_Base = MeanTime_Quadratic;
    using MeanTime_Hyporheic_Adv_Base = MeanTime_Linear;
    using HyporeicTransitionRegion_Base = HyporheicTransitionRegion_Probability;
    using HyporeicTransitionRegion_Adv_Base = HyporheicTransitionRegion_Probability;
    
    using MeanTime_WaterColumn_Surge = MeanTime_WaterColumn_Base;
    using MeanTime_Hyporheic_Adv_Surge = MeanTime_Hyporheic_Adv_Base;
    using HyporeicTransitionRegion_Surge = HyporeicTransitionRegion_Base;
    using HyporeicTransitionRegion_Adv_Surge = HyporheicTransitionRegion_Resuspension;
    
    using Time_Bed = stochastic::RNG<stochastic::skewedlevystable_distribution<double>>;

    MeanTime_WaterColumn_Base make_MeanTime_WaterColumn_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_w, 1./parameters.time_w_diff };
    }
    MeanTime_Hyporheic_Adv_Base make_MeanTime_Hyporheic_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_h  };
    }
    HyporeicTransitionRegion_Base make_HyporeicTransitionRegion_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp };
    }
    HyporeicTransitionRegion_Adv_Base make_HyporeicTransitionRegion_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp_adv };
    }
    
    MeanTime_WaterColumn_Surge make_MeanTime_WaterColumn_Surge
    (Parameters const& parameters)
    {
      return make_MeanTime_WaterColumn_Base(parameters);
    }
    MeanTime_Hyporheic_Adv_Surge make_MeanTime_Hyporheic_Adv_Surge
    (Parameters const& parameters)
    {
      return make_MeanTime_Hyporheic_Adv_Base(parameters);
    }
    HyporeicTransitionRegion_Surge make_HyporeicTransitionRegion_Surge
    (Parameters const& parameters)
    {
      return make_HyporeicTransitionRegion_Base(parameters);
    }
    HyporeicTransitionRegion_Adv_Surge make_HyporeicTransitionRegion_Adv_Surge
    (Parameters const& parameters)
    {
      return {};
    }

    Time_Bed make_Time_Bed
    (Parameters const& parameters)
    {
      return Time_Bed::param_type{
        parameters.alpha_b, parameters.time_b };
    }
  }
  
  namespace model_3
  {
    const std::string model_name = { "3" };
    
    struct Parameters
    {
      double time_w_diff;
      double time_h_diff;
      double time_b;
      double alpha_b;
      double advection_coeff_w;
      double advection_coeff_h;
      double prob_resusp;
      double prob_resusp_adv;
      double rate_decay_w;
      double rate_decay_h;
      double rate_decay_b;
      double concentration_inj_w;
      double mass_inj_h;
      double mass_inj_b;
      double time_injection_min;
      double time_injection_max;
      double time_step_injection;
      double max_dist;
      double time_min;
      double time_max;
      std::size_t nr_measurements;
      std::string measure_spacing;
      std::size_t nr_particles;
      std::string flow_name;
      std::size_t run_nr;
      std::string zones_name;
      
      Parameters(std::string const& filename)
      {
        auto input = useful::open_read(filename);
        useful::read(input, time_w_diff);
        useful::read(input, time_h_diff);
        useful::read(input, time_b);
        useful::read(input, alpha_b);
        useful::read(input, advection_coeff_w);
        useful::read(input, advection_coeff_h);
        useful::read(input, prob_resusp);
        useful::read(input, prob_resusp_adv);
        useful::read(input, rate_decay_w);
        useful::read(input, rate_decay_h);
        useful::read(input, rate_decay_b);
        useful::read(input, concentration_inj_w);
        useful::read(input, mass_inj_h);
        useful::read(input, mass_inj_b);
        useful::read(input, time_injection_min);
        useful::read(input, time_injection_max);
        useful::read(input, time_step_injection);
        useful::read(input, max_dist);
        useful::read(input, time_min);
        useful::read(input, time_max);
        useful::read(input, nr_measurements);
        useful::read(input, measure_spacing);
        useful::read(input, nr_particles);
        input.close();
      }
      
      static void parameter_list()
      {
        std::cout << "Model name: " << model_name << "\n"
                  << "Parameters (default in []):\n"
                  << "\ttime_w_diff             : Mean dispersive residence time in water column\n"
                  << "\ttime_h_diff             : Mean dispersive residence time in hyporheic zone\n"
                  << "\ttime_b                  : Characteristic residence time in bed\n"
                  << "\talpha_b                 : Residence time exponent in bed\n"
                  << "\tadvection_coeff_w       : Coefficient of advection-induced rate in water column\n"
                  << "\tadvection_coeff_h       : Coefficient of advection-induced rate in hyporheic zone\n"
                  << "\tprob_resusp             : Probability of transition to water column by diffusion during base flow\n"
                  << "\tprob_resusp_adv         : Probability of transition to water column by advection during base flow\n"
                  << "\trate_decay_w            : Decay rate in water column\n"
                  << "\trate_decay_h            : Decay rate in hyporheic zone\n"
                  << "\trate_decay_b            : Decay rate in bed\n"
                  << "\tconcentration_inj_w     : Concentration (mass if pulse) at upstream boundary in water column\n"
                  << "\tmass_inj_h              : Mass in pulse at upstream boundary in hyporheic zone\n"
                  << "\tmass_inj_b              : Mass in pulse at upstream boundary in bed\n"
                  << "\ttime_injection_min      : Start time of injection\n"
                  << "\ttime_injection_max      : End time of injection (equal to start time for pulse)\n"
                  << "\ttime_step_injection     : Discretization time step for injection (ignored for pulse)\n"
                  << "\tmax_dist                : Maximum downstream distance\n"
                  << "\ttime_min                : Minimum measure time\n"
                  << "\ttime_max                : Maximum measure time\n"
                  << "\tnr_measurements         : Nr of measurement times\n"
                  << "\tmeasure_spacing         : lin - Linear time spacing\n"
                  << "\t                          log - Logarithmic time spacing\n"
                  << "\tnr_particles            : Number of particles per injection\n";
      }
    };
    
    using MeanTime_WaterColumn_Base = MeanTime_Quadratic;
    using MeanTime_Hyporheic_Adv_Base = MeanTime_Linear;
    using HyporeicTransitionRegion_Base = HyporheicTransitionRegion_Probability;
    using HyporeicTransitionRegion_Adv_Base = HyporheicTransitionRegion_Probability;
    
    using MeanTime_WaterColumn_Surge = MeanTime_Infinite;
    using MeanTime_Hyporheic_Adv_Surge = MeanTime_Hyporheic_Adv_Base;
    using HyporeicTransitionRegion_Surge = HyporeicTransitionRegion_Base;
    using HyporeicTransitionRegion_Adv_Surge = HyporheicTransitionRegion_Resuspension;
    
    using Time_Bed = stochastic::RNG<stochastic::skewedlevystable_distribution<double>>;

    MeanTime_WaterColumn_Base make_MeanTime_WaterColumn_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_w, 1./parameters.time_w_diff };
    }
    MeanTime_Hyporheic_Adv_Base make_MeanTime_Hyporheic_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_h };
    }
    HyporeicTransitionRegion_Base make_HyporeicTransitionRegion_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp };
    }
    HyporeicTransitionRegion_Adv_Base make_HyporeicTransitionRegion_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp_adv };
    }
    
    MeanTime_WaterColumn_Surge make_MeanTime_WaterColumn_Surge
    (Parameters const& parameters)
    {
      return {};
    }
    MeanTime_Hyporheic_Adv_Surge make_MeanTime_Hyporheic_Adv_Surge
    (Parameters const& parameters)
    {
      return make_MeanTime_Hyporheic_Adv_Base(parameters);
    }
    HyporeicTransitionRegion_Surge make_HyporeicTransitionRegion_Surge
    (Parameters const& parameters)
    {
      return make_HyporeicTransitionRegion_Base(parameters);
    }
    HyporeicTransitionRegion_Adv_Surge make_HyporeicTransitionRegion_Adv_Surge
    (Parameters const& parameters)
    {
      return {};
    }

    Time_Bed make_Time_Bed
    (Parameters const& parameters)
    {
      return Time_Bed::param_type{
        parameters.alpha_b, parameters.time_b };
    }
  }
  
  namespace model_4
  {
    const std::string model_name = { "4" };
    
    struct Parameters
    {
      double time_w_diff;
      double time_h_diff;
      double time_b;
      double alpha_b;
      double advection_coeff_w;
      double advection_coeff_h;
      double advection_cutoff_h;
      double rate_cutoff_h;
      double prob_resusp;
      double prob_resusp_adv;
      double rate_decay_w;
      double rate_decay_h;
      double rate_decay_b;
      double concentration_inj_w;
      double mass_inj_h;
      double mass_inj_b;
      double time_injection_min;
      double time_injection_max;
      double time_step_injection;
      double max_dist;
      double time_min;
      double time_max;
      std::size_t nr_measurements;
      std::string measure_spacing;
      std::size_t nr_particles;
      
      Parameters(std::string const& filename)
      {
        auto input = useful::open_read(filename);
        useful::read(input, time_w_diff);
        useful::read(input, time_h_diff);
        useful::read(input, time_b);
        useful::read(input, alpha_b);
        useful::read(input, advection_coeff_w);
        useful::read(input, advection_coeff_h);
        useful::read(input, advection_cutoff_h);
        useful::read(input, rate_cutoff_h);
        useful::read(input, prob_resusp);
        useful::read(input, prob_resusp_adv);
        useful::read(input, rate_decay_w);
        useful::read(input, rate_decay_h);
        useful::read(input, rate_decay_b);
        useful::read(input, concentration_inj_w);
        useful::read(input, mass_inj_h);
        useful::read(input, mass_inj_b);
        useful::read(input, time_injection_min);
        useful::read(input, time_injection_max);
        useful::read(input, time_step_injection);
        useful::read(input, max_dist);
        useful::read(input, time_min);
        useful::read(input, time_max);
        useful::read(input, nr_measurements);
        useful::read(input, measure_spacing);
        useful::read(input, nr_particles);
        input.close();
      }
      
      static void parameter_list()
      {
        std::cout << "Model name: " << model_name << "\n"
                  << "Parameters (default in []):\n"
                  << "\ttime_w_diff             : Mean dispersive residence time in water column\n"
                  << "\ttime_h_diff             : Mean dispersive residence time in hyporheic zone\n"
                  << "\ttime_b                  : Characteristic residence time in bed\n"
                  << "\talpha_b                 : Residence time exponent in bed\n"
                  << "\tadvection_coeff_w       : Coefficient of advection-induced rate in water column\n"
                  << "\tadvection_coeff_h       : Coefficient of advection-induced rate in hyporheic zone\n"
                  << "\tadvection_cutoff_h      : Minimum value of advection for nonzero advection-induced rate in hyporheic zone\n"
                  << "\trate_cutoff_h           : Minimum value of advection-induced rate in hyporheic zone when advection exceeds minimum\n"
                  << "\tprob_resusp             : Probability of transition to water column by diffusion during base flow\n"
                  << "\tprob_resusp_adv         : Probability of transition to water column by advection during base flow\n"
                  << "\trate_decay_w            : Decay rate in water column\n"
                  << "\trate_decay_h            : Decay rate in hyporheic zone\n"
                  << "\trate_decay_b            : Decay rate in bed\n"
                  << "\tconcentration_inj_w     : Concentration (mass if pulse) at upstream boundary in water column\n"
                  << "\tmass_inj_h              : Mass in pulse at upstream boundary in hyporheic zone\n"
                  << "\tmass_inj_b              : Mass in pulse at upstream boundary in bed\n"
                  << "\ttime_injection_min      : Start time of injection\n"
                  << "\ttime_injection_max      : End time of injection (equal to start time for pulse)\n"
                  << "\ttime_step_injection     : Discretization time step for injection (ignored for pulse)\n"
                  << "\tmax_dist                : Maximum downstream distance\n"
                  << "\ttime_min                : Minimum measure time\n"
                  << "\ttime_max                : Maximum measure time\n"
                  << "\tnr_measurements         : Nr of measurement times\n"
                  << "\tmeasure_spacing         : lin - Linear time spacing\n"
                  << "\t                          log - Logarithmic time spacing\n"
                  << "\tnr_particles            : Number of particles per injection\n";
      }
    };
    
    using MeanTime_WaterColumn_Base = MeanTime_Quadratic;
    using MeanTime_Hyporheic_Adv_Base = MeanTime_Linear_Cutoff;
    using HyporeicTransitionRegion_Base = HyporheicTransitionRegion_Probability;
    using HyporeicTransitionRegion_Adv_Base = HyporheicTransitionRegion_Probability;
    
    using MeanTime_WaterColumn_Surge = MeanTime_WaterColumn_Base;
    using MeanTime_Hyporheic_Adv_Surge = MeanTime_Hyporheic_Adv_Base;
    using HyporeicTransitionRegion_Surge = HyporeicTransitionRegion_Base;
    using HyporeicTransitionRegion_Adv_Surge = HyporheicTransitionRegion_Resuspension;
    
    using Time_Bed = stochastic::RNG<stochastic::skewedlevystable_distribution<double>>;

    MeanTime_WaterColumn_Base make_MeanTime_WaterColumn_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_w, 1./parameters.time_w_diff };
    }
    MeanTime_Hyporheic_Adv_Base make_MeanTime_Hyporheic_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_h, parameters.advection_cutoff_h,
        parameters.rate_cutoff_h };
    }
    HyporeicTransitionRegion_Base make_HyporeicTransitionRegion_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp };
    }
    HyporeicTransitionRegion_Adv_Base make_HyporeicTransitionRegion_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp_adv };
    }
    
    MeanTime_WaterColumn_Surge make_MeanTime_WaterColumn_Surge
    (Parameters const& parameters)
    {
      return make_MeanTime_WaterColumn_Base(parameters);
    }
    MeanTime_Hyporheic_Adv_Surge make_MeanTime_Hyporheic_Adv_Surge
    (Parameters const& parameters)
    {
      return make_MeanTime_Hyporheic_Adv_Base(parameters);
    }
    HyporeicTransitionRegion_Surge make_HyporeicTransitionRegion_Surge
    (Parameters const& parameters)
    {
      return make_HyporeicTransitionRegion_Base(parameters);
    }
    HyporeicTransitionRegion_Adv_Surge make_HyporeicTransitionRegion_Adv_Surge
    (Parameters const& parameters)
    {
      return {};
    }

    Time_Bed make_Time_Bed
    (Parameters const& parameters)
    {
      return Time_Bed::param_type{
        parameters.alpha_b, parameters.time_b };
    }
  }
  
  namespace model_5
  {
    const std::string model_name = { "5" };
    
    struct Parameters
    {
      double time_w_diff;
      double time_h_diff;
      double time_b;
      double alpha_b;
      double advection_coeff_w;
      double advection_coeff_h;
      double prob_resusp;
      double prob_resusp_adv;
      double rate_decay_w;
      double rate_decay_h;
      double rate_decay_b;
      double concentration_inj_w;
      double mass_inj_h;
      double mass_inj_b;
      double time_injection_min;
      double time_injection_max;
      double time_step_injection;
      double max_dist;
      double time_min;
      double time_max;
      std::size_t nr_measurements;
      std::string measure_spacing;
      std::size_t nr_particles;
      
      Parameters(std::string const& filename)
      {
        auto input = useful::open_read(filename);
        useful::read(input, time_w_diff);
        useful::read(input, time_h_diff);
        useful::read(input, time_b);
        useful::read(input, alpha_b);
        useful::read(input, advection_coeff_w);
        useful::read(input, advection_coeff_h);
        useful::read(input, prob_resusp);
        useful::read(input, prob_resusp_adv);
        useful::read(input, rate_decay_w);
        useful::read(input, rate_decay_h);
        useful::read(input, rate_decay_b);
        useful::read(input, concentration_inj_w);
        useful::read(input, mass_inj_h);
        useful::read(input, mass_inj_b);
        useful::read(input, time_injection_min);
        useful::read(input, time_injection_max);
        useful::read(input, time_step_injection);
        useful::read(input, max_dist);
        useful::read(input, time_min);
        useful::read(input, time_max);
        useful::read(input, nr_measurements);
        useful::read(input, measure_spacing);
        useful::read(input, nr_particles);
        input.close();
      }
      
      static void parameter_list()
      {
        std::cout << "Model name: " << model_name << "\n"
                  << "Parameters (default in []):\n"
                  << "\ttime_w_diff             : Mean dispersive residence time in water column\n"
                  << "\ttime_h_diff             : Mean dispersive residence time in hyporheic zone\n"
                  << "\ttime_b                  : Characteristic residence time in bed\n"
                  << "\talpha_b                 : Residence time exponent in bed\n"
                  << "\tadvection_coeff_w       : Coefficient of advection-induced rate in water column\n"
                  << "\tadvection_coeff_h       : Coefficient of advection-induced rate in hyporheic zone\n"
                  << "\tprob_resusp             : Probability of transition to water column by diffusion during base flow\n"
                  << "\tprob_resusp_adv         : Probability of transition to water column by advection during base flow\n"
                  << "\trate_decay_w            : Decay rate in water column\n"
                  << "\trate_decay_h            : Decay rate in hyporheic zone\n"
                  << "\trate_decay_b            : Decay rate in bed\n"
                  << "\tconcentration_inj_w     : Concentration (mass if pulse) at upstream boundary in water column\n"
                  << "\tmass_inj_h              : Mass in pulse at upstream boundary in hyporheic zone\n"
                  << "\tmass_inj_b              : Mass in pulse at upstream boundary in bed\n"
                  << "\ttime_injection_min      : Start time of injection\n"
                  << "\ttime_injection_max      : End time of injection (equal to start time for pulse)\n"
                  << "\ttime_step_injection     : Discretization time step for injection (ignored for pulse)\n"
                  << "\tmax_dist                : Maximum downstream distance\n"
                  << "\ttime_min                : Minimum measure time\n"
                  << "\ttime_max                : Maximum measure time\n"
                  << "\tnr_measurements         : Nr of measurement times\n"
                  << "\tmeasure_spacing         : lin - Linear time spacing\n"
                  << "\t                          log - Logarithmic time spacing\n"
                  << "\tnr_particles            : Number of particles per injection\n";
      }
    };
    
    using MeanTime_WaterColumn_Base = MeanTime_Quadratic;
    using MeanTime_Hyporheic_Adv_Base = MeanTime_Quadratic;
    using HyporeicTransitionRegion_Base = HyporheicTransitionRegion_Probability;
    using HyporeicTransitionRegion_Adv_Base = HyporheicTransitionRegion_Probability;
    
    using MeanTime_WaterColumn_Surge = MeanTime_WaterColumn_Base;
    using MeanTime_Hyporheic_Adv_Surge = MeanTime_Hyporheic_Adv_Base;
    using HyporeicTransitionRegion_Surge = HyporeicTransitionRegion_Base;
    using HyporeicTransitionRegion_Adv_Surge = HyporheicTransitionRegion_Resuspension;
    
    using Time_Bed = stochastic::RNG<stochastic::skewedlevystable_distribution<double>>;

    MeanTime_WaterColumn_Base make_MeanTime_WaterColumn_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_w, 1./parameters.time_w_diff };
    }
    MeanTime_Hyporheic_Adv_Base make_MeanTime_Hyporheic_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.advection_coeff_h };
    }
    HyporeicTransitionRegion_Base make_HyporeicTransitionRegion_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp };
    }
    HyporeicTransitionRegion_Adv_Base make_HyporeicTransitionRegion_Adv_Base
    (Parameters const& parameters)
    {
      return { parameters.prob_resusp_adv };
    }
    
    MeanTime_WaterColumn_Surge make_MeanTime_WaterColumn_Surge
    (Parameters const& parameters)
    {
      return make_MeanTime_WaterColumn_Base(parameters);
    }
    MeanTime_Hyporheic_Adv_Surge make_MeanTime_Hyporheic_Adv_Surge
    (Parameters const& parameters)
    {
      return make_MeanTime_Hyporheic_Adv_Base(parameters);
    }
    HyporeicTransitionRegion_Surge make_HyporeicTransitionRegion_Surge
    (Parameters const& parameters)
    {
      return make_HyporeicTransitionRegion_Base(parameters);
    }
    HyporeicTransitionRegion_Adv_Surge make_HyporeicTransitionRegion_Adv_Surge
    (Parameters const& parameters)
    {
      return {};
    }

    Time_Bed make_Time_Bed
    (Parameters const& parameters)
    {
      return Time_Bed::param_type{
        parameters.alpha_b, parameters.time_b };
    }
  }
}

#endif /* Models_EColi_h */

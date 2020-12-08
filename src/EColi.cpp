//
// main.cpp
// EColi
//
// Created by Tomas Aquino on 8/2/19.
// Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Simulate E. coli dynamics in stream reach.
// E. coli moves by piecewise-constant-in-time advection in water column,
// no downstream movement elsewhere.
// Exponential residence time in hyporheic zone due to diffusion/advection,
// with diffusion leading to equal probability of transition to water column or bed.
// and advection leading to transition to water column only.
// Stable residence time bed.
// Linear decay with region-dependent rates
// Constant injection boundary condition upstream
// Absorbing boundary condition downstream

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "EColi/EColi.h"
#include "EColi/Measurer_EColi.h"
#include "EColi/Models.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "Stochastic/CTRW/CTRW.h"

int main(int argc, const char * argv[])
{
  using namespace ecoli::model_3;
  
  if (argc == 1)
  {
    Parameters::parameter_list();
    return 0;
  }

  Parameters parameters(argc, argv);
  
//  std::string flow_dir = "../flow";
//  std::string output_dir = "../data";
  
  std::string flow_dir = "/Users/tomasaquino/Dropbox/EColi/code/flow";
  std::string output_dir = "/Users/tomasaquino/Dropbox/EColi/code/data";

  auto time_watercolumn = make_MeanTime_WaterColumn_Base(parameters);
  auto time_hyporheic_adv = make_MeanTime_Hyporheic_Adv_Base(parameters);
  auto hyporheic_transition_region = make_HyporeicTransitionRegion_Base(parameters);
  auto hyporheic_transition_region_adv = make_HyporeicTransitionRegion_Adv_Base(parameters);
  auto time_watercolumn_surge = make_MeanTime_WaterColumn_Surge(parameters);
  auto time_hyporheic_adv_surge = make_MeanTime_Hyporheic_Adv_Surge(parameters);
  auto hyporheic_transition_region_surge = make_HyporeicTransitionRegion_Surge(parameters);
  auto hyporheic_transition_region_adv_surge = make_HyporeicTransitionRegion_Adv_Surge(parameters);
  auto time_bed = make_Time_Bed(parameters);

  using State = ecoli::State_EColi;
  using CTRW = ctrw::CTRW<State>;
  using Measurer = ecoli::Measurer_Store_Mass_Region;

  ecoli::Transitions_EColi transitions{
    { flow_dir + "/flow_" + parameters.flow_name + ".dat",
      flow_dir + "/surge_" + parameters.flow_name + ".dat" },  // Advection info
    { parameters.rate_decay_w,
      parameters.rate_decay_h, parameters.rate_decay_b },      // Reaction parameters
    parameters.max_dist,                                       // Distance to absorbing boundary
    ecoli::TransitionHelper_EColi_PiecewiseAdvection{
      parameters.time_h_diff,
      time_watercolumn, time_hyporheic_adv,
      hyporheic_transition_region, hyporheic_transition_region_adv,
      time_watercolumn_surge, time_hyporheic_adv_surge,
      hyporheic_transition_region_surge, hyporheic_transition_region_adv_surge,
      time_bed }
  };

  auto measure_times =
  [&parameters]()
  {
    if (parameters.measure_spacing.compare("lin") == 0)
      return range::linspace(parameters.time_min,
                             parameters.time_max, parameters.nr_measurements);
    if (parameters.measure_spacing.compare("log") == 0)
      return range::logspace(parameters.time_min,
                             parameters.time_max, parameters.nr_measurements);

    throw useful::bad_parameters();
  };
  
  // Prepare output
  std::string filename_output{ output_dir + "/Data_EColi_RegionMass_" + model_name + "_"
    + parameters.parameter_str() + ".dat" };
  
  std::cout << std::setprecision(2)
            << std::scientific;

  Measurer measurer{ measure_times() };
  
  // Handle continuous injection by discretizing
  // and simulating separately to conserve memory
  bool instantaneous = parameters.time_injection_max == parameters.time_injection_min
  ? 1 : 0;
  std::vector<double> injection_times = range::range(parameters.time_injection_min,
    parameters.time_step_injection, parameters.time_injection_max);
  for (auto time_injection : injection_times)
  {
    static bool first_injection = 1;
    double velocity = transitions.advection().velocity(time_injection);
    auto make_particles =
    [&parameters, velocity, instantaneous, time_injection]()
    {
      double factor = instantaneous ? 1. : velocity*parameters.time_step_injection;
      std::vector<double> mass{ factor*parameters.concentration_inj_w, 0., 0. };
      if (first_injection)
      {
        mass[1] = parameters.mass_inj_h;
        mass[2] = parameters.mass_inj_b;
        first_injection = 0;
      }
      double total_mass = operation::sum(mass);
      std::vector<std::size_t> particles_region = { 0,
        std::size_t(mass[1]/total_mass*parameters.nr_particles),
        std::size_t(mass[2]/total_mass*parameters.nr_particles) };
      particles_region[0] = parameters.nr_particles -
        (particles_region[1]+particles_region[2]);
      
      std::vector<CTRW::Particle> particles;
      particles.reserve(parameters.nr_particles);
      for (std::size_t region = 0; region < particles_region.size(); ++region)
        for (std::size_t pp = 0; pp < particles_region[region]; ++pp)
          particles.push_back(State{ 0.,
            mass[region]/particles_region[region],
            region, time_injection });
      
      return particles;
    };
    CTRW ctrw{ make_particles() };

    std::cout << "injection time = " << time_injection << "\n";
    for (auto time_measure : measurer.measure_times)
    {
      // Ignore measure times before current injection time
      if (time_measure < time_injection)
      {
        measurer.skip();
        continue;
      }

      // Evolve each particle to current measure time or until absorbed
      ctrw.evolve
      ([time_measure](auto const& part)
       { return part.state_new().time < time_measure
         && part.state_new().region != 3; },
       transitions);
      measurer(ctrw, transitions.reaction());
    }
    measurer.reset_measure();
  }
  measurer.print(filename_output);
  
  return 0;
}

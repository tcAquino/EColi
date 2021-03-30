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
  using namespace ecoli::model_1;
  
  if (argc == 1)
  {
    std::cout << "EColi : usage\n"
              << "--------------------------------------------------\n"
              << "Executable parameters (default in []):\n"
              << "\tparameters_name   : Name of parameter set. Read parameters\n"
              << "\t                    from parameters_parameters_name.dat\n"
              << "\tflow_name         : Name of flow data. Read flow data\n"
              << "\t                    from ../flow/flow_flow_name.dat,\n"
              << "\t                    and surge data (if applicable)\n"
              << "\t                    from ../flow/surge_flow_name.dat\n"
              << "\trun_nr            : Nonnegative integer index of output data file\n"
              << "\tzones_name [\"\"] : Name of measurement zone specifications. Read zone\n"
              << "\t                    specifications from ../data/zones_zones_name.dat\n"
              << "--------------------------------------------------\n"
              << "Remaining parameters, see below,\n"
              << "to be placed in parameters_parameters_name.dat\n"
              << "--------------------------------------------------\n";
    Parameters::parameter_list();
    return 0;
  }

  if (argc != 4 && argc != 5)
    throw useful::bad_parameters();
  
  std::size_t arg = 1;
  std::string parameters_name = argv[arg++];
  std::string flow_name = argv[arg++];
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string zones_name = argc > arg
  ? argv[arg++]
  : "";
  
  std::string parameters_dir = "../parameters";
  Parameters parameters(parameters_dir + "/parameters_" + parameters_name);
  
  std::string flow_dir = "../flow";
  std::string output_dir = "../data";

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
  using Measurer_Domain = ecoli::Measurer_Store_Mass_Region;
  using Measurer_Zones = ecoli::Measurer_Store_Mass_Region_Zones;

  ecoli::Transitions_EColi transitions{
    { flow_dir + "/flow_" + flow_name + ".dat",
      flow_dir + "/surge_" + flow_name + ".dat" },  // Advection info
    { parameters.rate_decay_w,
      parameters.rate_decay_h, parameters.rate_decay_b },           // Reaction parameters
    parameters.max_dist,                                            // Distance to absorbing boundary
    ecoli::TransitionHelper_EColi_PiecewiseAdvection{
      parameters.time_h_diff,
      time_watercolumn, time_hyporheic_adv,
      hyporheic_transition_region, hyporheic_transition_region_adv,
      time_watercolumn_surge, time_hyporheic_adv_surge,
      hyporheic_transition_region_surge, hyporheic_transition_region_adv_surge,
      time_bed }
  };

  auto get_measure_times =
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
  
  auto get_measure_zones =
  [&zones_name]()
  {
    std::vector<std::pair<double, double>> zones;
    if (zones_name.empty())
      return zones;
    
    auto file = useful::open_read("../data/zones_"
      + zones_name + ".dat");
    double left, right;
    while (file >> left >> right)
      zones.push_back({ left, right });
    return zones;
  };
  
  // Prepare output
  std::string filename_output_domain{ output_dir + "/"
    + "Data_EColi_RegionMass_" +
    + "model_" + model_name + "_"
    + "parameters_" + parameters_name + "_"
    + "flow_" + flow_name + "_"
    + "run_" + std::to_string(run_nr) + ".dat" };
  std::string filename_output_zones{ output_dir + "/"
    + "Data_EColi_RegionMass_Zones" +
    + "model_" + model_name + "_"
    + "parameters_" + parameters_name + "_"
    + "flow_" + flow_name + "_"
    + "zones_" + zones_name + "_"
    + "run_" + std::to_string(run_nr) + ".dat" };
  
  std::cout << std::setprecision(2)
            << std::scientific;

  auto measure_times = get_measure_times();
  Measurer_Domain measurer_domain{ measure_times };
  Measurer_Zones measurer_zones{ measure_times, get_measure_zones() };
  std::cout << measurer_zones.zones.size() << "\n";
  
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
    for (auto time_measure : measure_times)
    {
      // Ignore measure times before current injection time
      if (time_measure < time_injection)
      {
        measurer_domain.skip();
        measurer_zones.skip();
        continue;
      }

      // Evolve each particle to current measure time or until absorbed
      ctrw.evolve
      ([time_measure](auto const& part)
       { return part.state_new().time < time_measure
         && part.state_new().region != 3; },
       transitions);
      measurer_domain(ctrw, transitions.reaction());
      measurer_zones(ctrw, transitions.reaction(), transitions.advection());
    }
    measurer_domain.reset_measure();
    measurer_zones.reset_measure();
  }
  measurer_domain.print(filename_output_domain);
  measurer_zones.print(filename_output_zones);
  
  return 0;
}

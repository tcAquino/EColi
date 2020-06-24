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
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include "EColi/EColi.h"
#include "EColi/Measurer_EColi.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "general/useful.h"
#include "Stochastic/CTRW/CTRW.h"

int main(int argc, const char * argv[])
{
  if (argc == 1)
  {
    std::cout << "EColi : parameters:\n"
              << "\ttime_w_disp             : Mean dispersive residence time in water column\n"
              << "\ttime_h_disp             : Mean dispersive residence time in hyporheic zone\n"
              << "\ttime_b                  : Characteristic residence time in bed\n"
              << "\talpha_b                 : Residence time exponent in bed\n"
              << "\tadvection_coeff_w       : Coefficient of advection-induced rate in water column\n"
              << "\tadvection_coeff_h       : Coefficient of advection-induced rate in hyporheic zone\n"
              << "\tprob_resusp             : Probability of transition to water column by dispersion\n"
              << "\tprob_resusp_adv         : Probability of transition to water column by advection\n"
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
              << "\ttime_min                : Minimum measure timeDecay rate in water column\n"
              << "\ttime_max                : Maximum measure timeDecay rate in hyporheic zone\n"
              << "\tnr_measurements         : Nr of measurement times\n"
              << "\tmeasure_spacing         : lin - Linear time spacing\n"
              << "\t                          log - Logarithmic time spacing\n"
              << "\tnr_particles            : Number of particles per injection\n"
              << "\tflow_name               : Name of flow data\n"
              << "\trun_nr                  : Nonnegative integer index of output data file\n";
    return 0;
  }
  
  if (argc != 26)
    throw useful::bad_parameters();

  std::size_t arg = 1;
  double time_w_disp = atof(argv[arg++]);
  double time_h_disp = atof(argv[arg++]);
  double time_b = atof(argv[arg++]);
  double alpha_b = atof(argv[arg++]);
  double advection_coeff_w = atof(argv[arg++]);
  double advection_coeff_h = atof(argv[arg++]);
  double prob_resusp = atof(argv[arg++]);
  double prob_resusp_adv = atof(argv[arg++]);
  double rate_decay_w = atof(argv[arg++]);
  double rate_decay_h = atof(argv[arg++]);
  double rate_decay_b = atof(argv[arg++]);
  double concentration_inj_w = atof(argv[arg++]);
  double mass_inj_h = atof(argv[arg++]);
  double mass_inj_b = atof(argv[arg++]);
  double time_injection_min = atof(argv[arg++]);
  double time_injection_max = atof(argv[arg++]);
  double time_step_injection = atof(argv[arg++]);
  double max_dist = atof(argv[arg++]);
  double time_min = atof(argv[arg++]);
  double time_max = atof(argv[arg++]);
  std::size_t nr_measurements = strtoul(argv[arg++], NULL, 0);
  std::string measure_spacing = argv[arg++];
  std::size_t nr_particles = strtoul(argv[arg++], NULL, 0);
  std::string flow_name = argv[arg++];
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  
  std::string flow_dir = "../flow";
  std::string output_dir = "../data";
  
//  std::string flow_dir = "/Users/tomasaquino/Dropbox/EColi/code/flow"
//  std::string output_dir = "/Users/tomasaquino/Dropbox/EColi/code/data";

  // Full mean residence time in water column (deposition time)
  auto time_watercolumn = [time_w_disp, advection_coeff_w]
  (double advection, bool surge)
  {
    // If during surge event, no deposition
    if (surge)
      return std::numeric_limits<double>::infinity();
    // If not during surge event, rate proportional to advection
    else
      return 1./(1./time_w_disp + advection_coeff_w*advection*advection);
  };
  
  // Mean resuspension time by advection in hyporheic zone
  auto time_hyporheic_adv = [advection_coeff_h]
  (double advection, bool surge)
  {
    // If during surge event, instantaneous resuspension
    if (surge)
      return 0.;
    // If not during surge event, rate proportional to advection
    else
      return 1./(advection_coeff_h*advection);
  };

  using State = ecoli::State_EColi;
  using CTRW = ctrw::CTRW<State>;
  using Time_WaterColumn = decltype(time_watercolumn);
  using Time_Hyporheic = decltype(time_hyporheic_adv);
  using TransitionHelper =
    ecoli::TransitionHelper_EColi_PiecewiseAdvection
    <Time_WaterColumn, Time_Hyporheic>;
  using Transitions = ecoli::Transitions_EColi<TransitionHelper>;
  using Measurer = ecoli::Measurer_Store_Mass_Region;

  Transitions transitions{
    { flow_dir + "/flow_" + flow_name + ".dat",
      flow_dir + "/surge_" + flow_name + ".dat"
    },                                               // Piecewise advection timeseries with surge info
    { time_h_disp, time_b, alpha_b,
      time_watercolumn, time_hyporheic_adv,
      prob_resusp, prob_resusp_adv },                // Residence time parameters
    { rate_decay_w, rate_decay_h, rate_decay_b },    // Reaction parameters
    max_dist };                                      // Distance to absorbing boundary

  auto measure_times =
  [measure_spacing, time_min, time_max, nr_measurements]()
  {
    if (measure_spacing.compare("lin") == 0)
      return range::linspace(time_min, time_max, nr_measurements);
    if (measure_spacing.compare("log") == 0)
      return range::logspace(time_min, time_max, nr_measurements);

    throw useful::bad_parameters();
  };
  
  // Prepare output
  std::cout << std::setprecision(2)
            << std::scientific;
  std::string filename_base{ output_dir + "/Data_EColi_RegionMass" };
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << time_w_disp << "_"
         << time_h_disp << "_"
         << time_b << "_"
         << alpha_b << "_"
         << advection_coeff_w << "_"
         << advection_coeff_h << "_"
         << prob_resusp << "_"
         << prob_resusp_adv << "_"
         << rate_decay_w << "_"
         << concentration_inj_w << "_"
         << mass_inj_h << "_"
         << mass_inj_b << "_"
         << time_injection_min << "_"
         << time_injection_max << "_"
         << time_step_injection << "_"
         << max_dist << "_"
         << time_min << "_"
         << time_max << "_"
         << nr_measurements << "_"
         << measure_spacing << "_"
         << nr_particles << "_"
         << flow_name << "_"
         << run_nr;
  std::string filename_params = stream.str();
  std::string filename_output = filename_base +
    filename_params + ".dat";

  Measurer measurer{ measure_times() };
  
  // Handle continuous injection by discretizing
  // and simulating separately to conserve memory
  bool instantaneous = time_injection_max == time_injection_min
  ? 1 : 0;
  std::vector<double> injection_times = range::range(time_injection_min,
    time_step_injection, time_injection_max);
  for (auto time_injection : injection_times)
  {
    static bool first_injection = 1;
    double velocity = transitions.advection().velocity(time_injection);
    auto make_particles =
    [nr_particles, time_injection, instantaneous,
     concentration_inj_w, mass_inj_h, mass_inj_b,
     velocity, time_step_injection]()
    {
      double factor = instantaneous ? 1. : velocity*time_step_injection;
      std::vector<double> mass{ factor*concentration_inj_w, 0., 0. };
      if (first_injection)
      {
        mass[1] = mass_inj_h;
        mass[2] = mass_inj_b;
        first_injection = 0;
      }
      double total_mass = operation::sum(mass);
      std::vector<std::size_t> particles_region = { 0,
        std::size_t(mass[1]/total_mass*nr_particles),
        std::size_t(mass[2]/total_mass*nr_particles) };
      particles_region[0] = nr_particles -
        (particles_region[1]+particles_region[2]);
      
      std::vector<CTRW::Particle> particles;
      particles.reserve(nr_particles);
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
      ([time_measure, max_dist](auto const& part)
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

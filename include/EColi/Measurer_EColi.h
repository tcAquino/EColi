//
//  Measurer_EColi_region.h
//  EColi
//
//  Created by Tomas Aquino on 8/5/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#ifndef Observer_region_h
#define Observer_region_h

#include <iomanip>
#include <valarray>
#include <vector>
#include "general/Operations.h"
#include "general/useful.h"

namespace ecoli
{
  // Measure the mass in each region
  // at each measure time
	class Measurer_Store_Mass_Region
	{
	public:
		Measurer_Store_Mass_Region
    (std::vector<double> measure_times)
    : measure_times{ measure_times }
    , masses(measure_times.size(), std::vector<double>(4))
    {}
      
    // Add values to the current measure time
    template <typename Subject, typename Reaction>
    void operator()
    (Subject const& subject, Reaction const& reaction)
    {
      for (auto const& part : subject.particles())
      {
        auto const& state_old = part.state_old();
        auto const& state_new = part.state_new();
        
        // Ignore particles with injection time which has not yet been reached
        if (state_old.time > measure_times[measure])
          continue;
          
        // If particle has not yet been absorbed in the new state
        if (state_new.region != 3
            && state_old.time < measure_times[measure])
        {
          interpolate(state_old, reaction);
          continue;
        }
        
        // If particle has been absorbed in the new state
        // If the new state is in the past
        if (state_new.time < measure_times[measure])
          masses[measure][3] += part.state_new().mass;
        // If it is in the future
        else
          interpolate(state_old, reaction);
      }
      ++measure;
    }
    
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      auto output = useful::open_write(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t tt = 0; tt < measure_times.size(); ++tt)
      {
        output << measure_times[tt];
        useful::print(output, masses[tt], 1, delimiter);
        output << "\n";
      }
      output.close();
    }
    
    // Skip current measure time without measuring anything
    void skip()
    { ++measure; }
    
    // Restart measuring from first measure time, keeping current values
    void reset_measure()
    { measure = 0; }
      
    const std::vector<double> measure_times;
    
	private:
    // masses[tt][rr] is the mass in region rr at time measure_times[tt]
    std::vector<std::vector<double>> masses;
    std::size_t measure{ 0 };
    
    // React on a copy of the old state until the measurement time
    template <typename State, typename Reaction>
    void interpolate
    (State const& state_old, Reaction const& reaction)
    {
      auto state = state_old;
      reaction(state, measure_times[measure]-state.time);
      masses[measure][state.region] += state.mass;
    }
	};
  
  // Measure the mass in each region within given downstream zones
  // at each measure time
  class Measurer_Store_Mass_Region_Zones
  {
  public:
    Measurer_Store_Mass_Region_Zones
    (std::vector<double> measure_times,
     std::vector<std::pair<double,double>> zones)
    : measure_times{ measure_times }
    , zones{ zones }
    , masses(measure_times.size(),
             std::vector<std::vector<double>>(
               zones.size(), std::vector<double>(3)))
    {}
      
    // Add values to the current measure time
    template <typename Subject, typename Reaction, typename Advection>
    void operator()
    (Subject const& subject, Reaction const& reaction,
     Advection const& advection)
    {
      for (auto const& part : subject.particles())
      {
        auto const& state_old = part.state_old();
        auto const& state_new = part.state_new();
        // Ignore particles with injection time which has not yet been reached
        if (state_old.time > measure_times[measure])
          continue;
          
        // If particle has not yet been absorbed in the new state
        if (state_new.region != 3
            && state_old.time < measure_times[measure])
        {
          interpolate(state_old, reaction, advection);
          continue;
        }
        
        // If particle has been absorbed in the new state
        // If the new state is in the future
        if (state_new.time > measure_times[measure])
          for (std::size_t zz = 0; zz < zones.size(); ++zz)
            if (state_new.position > zones[zz].first
                && state_new.position < zones[zz].second)
              interpolate(state_old, reaction, advection);
      }
      ++measure;
    }
    
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      if (zones.size() == 0)
        return;
      auto output = useful::open_write(filename);
      output << std::setprecision(precision)
             << std::scientific;
      for (std::size_t tt = 0; tt < measure_times.size(); ++tt)
      {
        output << measure_times[tt];
        for (auto const& mass_zone : masses[tt])
          useful::print(output, mass_zone, 1, delimiter);
        output << "\n";
      }
      output.close();
    }
    
    // Skip current measure time without measuring anything
    void skip()
    { ++measure; }
    
    // Restart measuring from first measure time, keeping current values
    void reset_measure()
    { measure = 0; }
      
    const std::vector<double> measure_times;
    const std::vector<std::pair<double,double>> zones;
    
  private:
    // masses[tt][zz][rr] is the mass in region rr within zone zz
    // at time measure_times[tt]
    std::vector<std::vector<std::vector<double>>> masses;
    std::size_t measure{ 0 };
    
    // React on and advect a copy of the old state until the measurement time
    template <typename State, typename Reaction, typename Advection>
    void interpolate
    (State const& state_old, Reaction const& reaction,
     Advection const& advection)
    {
      auto state = state_old;
      double exposure_time = measure_times[measure]-state.time;
      reaction(state, exposure_time);
      operation::plus_InPlace(state.position, advection.jump(state, exposure_time));
      for (std::size_t zz = 0; zz < zones.size(); ++zz)
        if (state.position >= zones[zz].first
            && state.position < zones[zz].second)
          masses[measure][zz][state.region] += state.mass;
    }
  };
}


#endif /* Observer_region_h */

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
    : measure_times(measure_times)
    , masses(measure_times.size(), std::vector<double>(4))
    {}
      
    // Each call adds values to the current measure time
    template <typename Subject, typename Reaction>
    void operator()
    (Subject const& subject, Reaction const& reaction)
    {
      for (auto const& part : subject.particles())
      {
        // Ignore particles with initial time which has not yet been reached
        if (part.state_old().time > measure_times[measure])
          continue;
        // If particle has not yet been aborbed
        if (part.state_new().region != 3)
        {
          // React on a copy of the old state until the measurement time
          auto state{ part.state_old() };
          reaction(state, measure_times[measure]-state.time);
          masses[measure][state.region] += state.mass;
        }
        // If particle has been absorbed
        else
          masses[measure][3] += part.state_new().mass;
      }
      ++measure;
    }
    
    void print
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t") const
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
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
    
    // Restart measuring from first measure time, keeping current values
    void reset_measure()
    { measure = 0; }
      
    const std::vector<double> measure_times;
    
	private:
    // masses[tt][rr] is the mass in region rr at time measure_times[tt]
    std::vector<std::vector<double>> masses;
    std::size_t measure{ 0 };
	};
}


#endif /* Observer_region_h */

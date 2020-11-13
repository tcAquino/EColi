//
//  DepositionResuspension.h
//  EColi
//
//  Created by Tomás Aquino on 01/11/2020.
//  Copyright © 2020 Tomas Aquino. All rights reserved.
//

#ifndef DepositionResuspension_h
#define DepositionResuspension_h

#include <limits>
#include <random>

namespace ecoli
{
  struct MeanTime_Quadratic
  {
    const double rate_base;
    const double advection_coeff;
    
    MeanTime_Quadratic(double rate_base, double advection_coeff)
    : rate_base{ rate_base }
    , advection_coeff{ advection_coeff }
    {}
    
    double operator() (double advection)
    {
      return 1./(rate_base + advection_coeff*advection*advection);
    }
  };

  struct MeanTime_Linear
  {
    const double advection_coeff;
    const double rate_base;
    
    MeanTime_Linear(double advection_coeff, double rate_base = 0.)
    : advection_coeff{ advection_coeff }
    , rate_base{ rate_base }
    {}
    
    double operator() (double advection)
    {
      return 1./(rate_base + advection_coeff*advection);
    }
  };

  struct MeanTime_LinearQuadratic
  {
    const double advection_coeff_linear;
    const double advection_coeff_quadratic;
    const double rate_base;
    
    MeanTime_LinearQuadratic
    (double advection_coeff_linear, double advection_coeff_quadratic,
     double rate_base = 0.)
    : advection_coeff_linear{ advection_coeff_linear }
    , advection_coeff_quadratic{ advection_coeff_quadratic }
    , rate_base{ rate_base }
    {}
    
    double operator() (double advection)
    {
      return 1./(rate_base + advection_coeff_linear*advection
                 + advection_coeff_quadratic*advection*advection);
    }
  };

  struct MeanTime_Linear_Cutoff
  {
    const double advection_coeff;
    const double min_val_adv;
    const double rate_base;
    
    MeanTime_Linear_Cutoff
    (double advection_coeff, double min_val_adv, double rate_base = 0.)
    : advection_coeff{ advection_coeff }
    , min_val_adv{ min_val_adv }
    , rate_base{ rate_base }
    {}
    
    double operator() (double advection)
    {
      double rate = advection > min_val_adv
      ? advection_coeff*(advection-min_val_adv)+rate_base
      : 0.;
      return 1./rate;
    }
  };

  struct MeanTime_Quadratic_Cutoff
  {
    const double advection_coeff;
    const double min_val_adv;
    const double rate_base;
    
    MeanTime_Quadratic_Cutoff
    (double advection_coeff, double min_val_adv, double rate_base = 0.)
    : advection_coeff{ advection_coeff }
    , min_val_adv{ min_val_adv }
    , rate_base{ rate_base }
    {}
    
    double operator() (double advection)
    {
      double rate = advection > min_val_adv
      ? advection_coeff*(advection*advection-min_val_adv*min_val_adv)+rate_base
      : 0.;
      return 1./rate;
    }
  };

  struct MeanTime_LinearQuadratic_Cutoff
  {
    const double advection_coeff_linear;
    const double advection_coeff_quadratic;
    const double min_val_adv;
    const double rate_base;
      
      MeanTime_LinearQuadratic_Cutoff
      (double advection_coeff_linear, double advection_coeff_quadratic,
      double min_val_adv, double rate_base = 0.)
      : advection_coeff_linear{ advection_coeff_linear }
      , advection_coeff_quadratic{ advection_coeff_quadratic }
      , min_val_adv{ min_val_adv }
      , rate_base{ rate_base }
      {}
      
      double operator() (double advection)
      {
        double rate_adv = advection > min_val_adv
        ? advection_coeff_linear*advection
          + advection_coeff_quadratic*advection*advection
        : 0.;
        return 1./(rate_base + rate_adv);
      }
  };
  
  struct MeanTime_Infinite
  {
    double operator() (double advection = 0.)
    { return std::numeric_limits<double>::infinity(); }
  };
  
  struct MeanTime_Zero
  {
    double operator() (double advection = 0.)
    { return 0.; }
  };
  
  struct MeanTime_Constant
  {
    const double mean_time;
    
    MeanTime_Constant(double mean_time)
    : mean_time{ mean_time }
    {}
    
    double operator() (double advection = 0.)
    { return mean_time; }
  };
  
  struct HyporheicTransitionRegion_Resuspension
  {
    std::size_t operator() ()
    { return 0; }
  };
  
  struct HyporheicTransitionRegion_Deposition
  {
    std::size_t operator() ()
    { return 2; }
  };
  
  struct HyporheicTransitionRegion_Probability
  {
    const double prob_resusp;
    
    HyporheicTransitionRegion_Probability(double prob_resusp)
    : prob_resusp{ prob_resusp }
    {}
    
    std::size_t operator() ()
    { return 2*dist_ber_diff(rng); }
    
  private:
    std::bernoulli_distribution dist_ber_diff{
      1.-prob_resusp };
    std::mt19937 rng{ std::random_device{}() };
  };
}


#endif /* DepositionResuspension_h */

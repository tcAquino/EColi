//
// Random.h
// Stochastic
//
// Created by Tomas Aquino on 2/7/15.
// Copyright (c) 2015 Tomas Aquino. All rights reserved.
//

// Random number generation and related utilities

#ifndef Stochastic_Random_h
#define Stochastic_Random_h

#include <array>
#include <cmath>
#include <random>
#include <utility>
#include <stdexcept>
#include <vector>
#include "general/useful.h"
#include "general/Constants.h"
#include "general/Operations.h"

namespace stochastic
{
  // Wrapper for random number generation with shared rng
  template <typename Distribution_t, typename Engine_t = std::mt19937>
  struct RNG_shared_engine
  {
    using param_type = typename Distribution_t::param_type;
    using result_type = typename Distribution_t::result_type;

    template <typename param_t>
    RNG_shared_engine(param_t const& params, Engine_t& rng)
    : dist{ param_type{ params } }
    , rng{ rng }
    {}

    result_type operator() ()
    { return dist(rng); }

    Distribution_t dist;

  private:
    Engine_t& rng;
  };

  // Wrapper for random number generation with own rng
  template <typename Distribution_t, typename Engine_t = std::mt19937>
  struct RNG
  {
    using param_type = typename Distribution_t::param_type;
    using result_type = typename Distribution_t::result_type;

    template <typename param_t>
    RNG(param_t const& params)
    : dist{ param_type(params) }
    {}

    result_type operator() ()
    { return dist(rng); }

    Distribution_t dist;

  private:
    Engine_t rng{ std::random_device{}() };
  };

  // Skewed Levy stable distribution
  // alpha : exponent, pdf $\sim t^{-1-\alpha}$
  // sigma : scale parameter
  // mu : location parameter
  // Note: tail \sim -sigma * tan(alpha*pi/2) * t^{-1-\alpha}/gamma(-alpha)
  template <typename Value_type = double>
  class skewedlevystable_distribution
  {
  public:
    using param_type = std::array<Value_type, 3>;
    using result_type = Value_type;

    const Value_type alpha;
    const Value_type sigma;
    const Value_type mu;

    const double zeta{ std::tan(constants::pi * alpha * 0.5) };
    const double xi{ std::atan(zeta) / alpha };
    const double vv{ std::pow(1. + zeta * zeta, 0.5 / alpha) };

    skewedlevystable_distribution(Value_type alpha, Value_type sigma = 1., Value_type mu = 0.)
    : alpha(alpha)
    , sigma(sigma)
    , mu(mu)
    {}

    skewedlevystable_distribution(param_type const& params)
    : alpha(params[0])
    , sigma(params[1])
    , mu(params[2])
    {}

    template <typename Generator>
    Value_type operator() (Generator& rng)
    {
      double uu = constants::pi * ( uniform_dist(rng) - 0.5 );
      double tt = std::sin(alpha * ( uu + xi ) ) / std::pow(std::cos( uu ), 1./alpha);
      double ss = std::pow(
                           std::cos(( 1. - alpha ) * uu - alpha * xi) / exponential_dist(rng),
                           ( 1.-alpha ) / alpha);

      return sigma * vv * tt * ss + mu;
    }

  private:
    std::uniform_real_distribution<double> uniform_dist{ 0., 1. };
    std::exponential_distribution<double> exponential_dist{ 1. };
  };

  template <typename Value_type = double>
  class pareto_distribution
  {
  public:
    using param_type = std::array<Value_type, 2>;
    using result_type = Value_type;

    const Value_type alpha;
    const Value_type min;

    pareto_distribution(Value_type alpha, Value_type min)
    : alpha(alpha)
    , min(min)
    {}

    pareto_distribution(param_type const& params)
    : alpha(params[0])
    , min(params[1])
    {}

    template <typename Generator>
    Value_type operator() (Generator& rng)
    {
      return min * std::exp(exponential_dist(rng) / alpha);
    }

  private:
    std::exponential_distribution<double> exponential_dist{ 1. };
  };

  // Inverse Gaussian distribution (with drift)
  // mean : dist mean
  // variance : dist variance
  template <typename Value_type = double>
  class inverse_gaussian_distribution
  {
  public:
    using param_type = std::array<Value_type, 2>;
    using result_type = Value_type;

    const Value_type mean;
    const Value_type var;

    inverse_gaussian_distribution(Value_type mean, Value_type var)
    : mean(mean)
    , var(var)
    {}

    inverse_gaussian_distribution(param_type const& params)
    : mean(params[0])
    , var(params[1])
    {}

    template <typename Generator>
    Value_type operator() (Generator& rng)
    {
      double rand_normal = normal_dist(rng);
      double aux = ratio * rand_normal * rand_normal;

      double xx = 1. + aux - std::sqrt(aux * ( aux + 2. ));
      double rand_uniform = uniform_dist(rng);

      return
      1. / ( 1. + xx ) < rand_uniform
      ? mean / xx
      : mean * xx;
    }

  private:
    std::uniform_real_distribution<double> uniform_dist{ 0., 1. };
    std::normal_distribution<double> normal_dist{ 0., 1. };

    double ratio{ 0.5*var/(mean*mean) };
  };
  
  template <typename Value_type = std::vector<double>>
  class isotropic_unit_vector_distribution
  {
  public:
    using param_type = std::size_t;
    using result_type = Value_type;

    isotropic_unit_vector_distribution(std::size_t dim)
    : dim{ dim }
    {}

    template <typename Generator>
    Value_type operator() (Generator& rng)
    {
      Value_type val;
      
      for(std::size_t dd = 0; dd < dim; ++dd)
        val.push_back(normal_dist(rng));
      operation::div_scalar_InPlace(val, operation::abs(val));
      return val;
    }

  private:
    param_type dim;
    std::normal_distribution<double> normal_dist{};
  };

  // Returns a pdf of given set of values,
  // normalized to integral equal to number of samples
  // Output bin values are bin centers
  // Data smaller than leftmost edge and larger than rightmost edge
  // are included in leftmost and rightmost bins, respectively
  template <typename Container>
  std::vector<std::pair<double, double>> pdf
  (std::vector<double> const& bin_edges, Container const& data)
  {
    std::vector<std::pair<double, double>> pdf(bin_edges.size() - 1);

    // Populate all bins except last
    auto rt_it = data.cbegin();
    for (std::size_t bin = 0; bin < bin_edges.size() - 2; ++bin)
    {
      pdf[bin].first = ( bin_edges[bin] + bin_edges[bin+1] ) / 2.;
      while (rt_it != data.cend()
             && *rt_it < bin_edges[bin+1])
      {
        ++pdf[bin].second;
        ++rt_it;
      }
      pdf[bin].second /= bin_edges[bin+1] - bin_edges[bin];
    }

    // Put everything that is left on last bin
    for (; rt_it != data.cend(); ++rt_it)
    {
      pdf.back().first = (bin_edges[bin_edges.size()-2] + bin_edges.back())/2.;
      ++pdf.back().second;
    }
    pdf.back().second /= bin_edges.back() - bin_edges[bin_edges.size()-2];

    return pdf;
  }

  // probs are cumulative probabilities (not necessarily normalized), p(0), p(0)+p(1), ....
  // Pick returns i with probability p(i)
  template<typename Container, typename Engine_t = std::mt19937>
  std::size_t pick(Container const& probs, Engine_t& rng)
  {
    typename Container::value_type rnd =
    probs.back() * std::uniform_real_distribution<double>{}(rng);

    for (std::size_t ev = 0; ev < probs.size(); ++ev)
      if (rnd < probs[ev]) return ev;

    throw std::runtime_error("Nothing picked!");
  }

  // Generates a random non-repeating sequence of
  // subset.size() elements out of {0, 1, ..., total-1}
  template<typename Engine_t = std::mt19937>
  void randSubset(std::size_t total, std::vector<std::size_t>& subset, Engine_t& rng)
  {
    std::vector<std::size_t> set;
    set.reserve(total);
    for (std::size_t ii = 0; ii < total; ++ii)
      set.push_back( ii );

    std::size_t picked;
    for (auto& ii : subset)
    {
      picked = std::uniform_int_distribution<std::size_t>{ 0, set.size()-1 }(rng);
      ii = set[picked];
      useful::swap_erase(set, picked);
    }
  }
  
  // Generates a random non-repeating sequence of
  // subset_size elements out of {0, 1, ..., total - 1}
  template<typename Engine_t = std::mt19937>
  std::vector<std::size_t> randSubset(std::size_t total, std::size_t subset_size, Engine_t& rng)
  {
    std::vector<std::size_t> subset(subset_size);
    RandSubset(total, subset, rng);
    
    return subset;
  }
}

#endif

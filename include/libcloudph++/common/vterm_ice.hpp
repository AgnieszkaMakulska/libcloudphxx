#pragma once

#include "units.hpp"
#include "macros.hpp"
#include "earth.hpp"
#include "const_cp.hpp"
#include "kelvin_term.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace vterm
    {

      // terminal fall velocity of ice based on Shima (2020), Eq. 5-6 and Eq. 111-121
      // formulas from Bohm (1992)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt_bohm(
        quantity<si::length, real_t> a, // equatorial radius
        quantity<si::length, real_t> c,  // polar radius
        quantity<si::mass_density, real_t> rho_i, // apparent density
        quantity<si::mass_density, real_t> rho, //density of moist air
        quantity<si::dynamic_viscosity, real_t> eta // air viscosity
      ) {

        assert(a > real_t(0) * si::meters);
        assert(c > real_t(0) * si::meters);

        using moist_air::rho_i_true;
        using earth::g;
#if !defined(__NVCC__)
        const real_t pi_ = pi<real_t>();
#else
        const real_t pi_ = CUDART_PI;
#endif

        real_t phi = c / a;
        quantity<si::mass, real_t> mi = real_t(4./3.) * pi_ * a * a * c * rho_i;

        // characteristic length
        quantity<si::length, real_t> d = real_t(2) * a;

        // area ratio
        quantity<si::dimensionless, real_t> q = std::pow(real_t(rho_i / rho_i_true), std::exp(- phi));

        // Davies number
        real_t X = real_t(8) * mi * g * rho /
          (pi * eta * eta * std::max(phi, real_t(1)) * std::pow(q, real_t(1./4.)));
        real_t X_ratio = X / real_t(2.8e6);
        real_t X_prime = X * ( real_t(1) + X_ratio * X_ratio)
          / (real_t(1) + real_t(1.6) * X_ratio * X_ratio);

        // viscous shape factor
        real_t k_1 = real_t(0.82) + real_t(0.18) * phi > real_t(0.85) ?
                     real_t(0.82) + real_t(0.18) * phi : real_t(0.85);
        real_t k_2 = real_t(0.37) + real_t(0.63) / std::sqrt(phi);
        real_t k_3 = real_t(1.33) / ((std::log(phi) > real_t(0) ? std::log(phi) : real_t(0)) + real_t(1.19));
        real_t k = k_1 < k_2 ? k_1 : k_2;
        k = k_3 < k ? k_3 : k;

        // gamma for the pressure drag coefficient
        real_t Gamma = real_t(3.76) - real_t(8.41) * phi + real_t(9.18) * phi * phi - real_t(3.53) * phi * phi * phi;
        Gamma = Gamma < real_t(1.98) ? Gamma : real_t(1.98);
        Gamma = Gamma > real_t(1) ? Gamma : real_t(1);

        // pressure drag coefficients
        real_t C_dp1 = real_t(0.292) * k * gamma;
        real_t C_dp2 = real_t(0.492) - real_t(0.2) / std::sqrt(phi);
        real_t C_dp = C_dp1 > C_dp2 ? C_dp1 : C_dp2;
        real_t C_do = real_t(4.5) * k * k * (phi > real_t(1) ? phi : real_t(1));

        real_t beta = std::sqrt(real_t(1) + C_dp / 6 / k * std::sqrt(X_prime / C_dp)) - real_t(1);
        real_t gam = (C_do - C_dp) / real_t(4) / C_dp;

        real_t N_re = real_t(6) * k * beta * beta / C_dp
          * (real_t(1) + real_t(2) * beta * std::exp(- beta * gam) / (real_t(2) + beta) / (real_t(1) + beta));

        return eta * N_re / (rho * d);
      }


    };
  };
};

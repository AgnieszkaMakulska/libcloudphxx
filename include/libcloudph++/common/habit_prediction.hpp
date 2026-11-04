#pragma once

#include "units.hpp"
#include "const_cp.hpp"
#include "detail/inherent_growth_ratio.hpp"

#include <libcloudph++/common/ventil.hpp>
#include <libcloudph++/common/vterm.hpp>

#include <boost/config.hpp>

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace common
  {
    namespace habit_prediction
    {

      // Depositional growth rate: dm/dt following Shima et al. 2020
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dmdt_i(
        const quantity<diffusivity, real_t> D,            // D
        const quantity<thermal_conductivity, real_t> K,   // K
        const quantity<si::mass_density, real_t> rho_v,   // ambient water vapour density
        const quantity<si::temperature, real_t> T,        // ambient temperature
        const quantity<si::dimensionless, real_t> RH_i,   // p_v/p_vsi
        const quantity<si::length, real_t> C,             // capacitance of the crystal
        const quantity<si::dimensionless, real_t> f_vnt   // ventilation coeff.
      )
      {
        using moist_air::R_v;

        quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t>
          l_s = const_cp::l_s<real_t>(T);

        return real_t(4) * pi<real_t>() * C * f_vnt
          * (RH_i - real_t(1))
          / (
            RH_i / D / rho_v                                   // F_d
            +
            l_s / K / T * (l_s / R_v<real_t>() / T - real_t(1)) // F_k
          );
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t capacitance(
        const real_t a,  // equatorial radius [m]
        const real_t c   // polar radius [m]
      )
      {
        const real_t phi =  c / a;
        if (phi >= real_t(1.001)) // prolate spheroid
        {
          const real_t ecent = std::sqrt(real_t(1.0) - real_t(1.0) / (std::pow(phi, 2)));
          return ecent * c / std::log((real_t(1) + ecent) * phi);
        }
        else if (phi <= real_t(0.999)) // oblate spheroid
        {
          const real_t ecent = std::sqrt(real_t(1) - std::pow(phi,2));
          return ecent * a / std::asin(ecent);
        }
        else // approximation around phi = 1
          return (real_t(2) * a + c) / real_t(3);
      }


      // Inherent growth ratio Gamma(T) of Chen and Lamb (1994)
      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t growth_ratio(
        const real_t T,  // temperature [K]
        const real_t a,  // equatorial radius [m]
        const real_t c   // polar radius [m]
      )
      {
        // quasi-spherical small crystals: Gamma = 1 if 2*max(a,c) < 10 um
        const real_t Dmax = real_t(2) * (a > c ? a : c);
        if (Dmax < real_t(1e-5)) return real_t(1);

        const real_t T_deg = T - real_t(273.15);

        // nearest table entry
        int n = int(real_t(1.375) - real_t(4) * T_deg + real_t(0.5));
        if (real_t(1.375) - real_t(4) * T_deg + real_t(0.5) < real_t(1)) n = 1;
        if (n < 1)   n = 1;
        if (n > 121) n = 121;

        return detail::inherent_growth_ratio_table<real_t>(n);
      }


      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t ice_ventilation_coeff(
      const real_t a,          // equatorial radius [m]
      const real_t c,          // polar radius [m]
      const real_t vterm_ice,  // terminal velocity of ice particle [m/s]
      const real_t D,          // vapor diffusivity
      const real_t eta,        // air viscosity
      const real_t rho         // moist air density [kg/m^3]
      ) {
        const real_t r_max = a > c ? a : c; // particle's max radius
        const real_t X = std::cbrt(ventil::Sc(eta, rho, D)) * std::sqrt(ventil::Re(vterm_ice, r_max, rho, eta));

        return (X <= real_t(1))
          ? real_t(1.0)  + real_t(0.14) * X * X
          : real_t(0.86) + real_t(0.28) * X;
      }


      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t ice_ventilation_coeff_ratio(
      const real_t a,          // equatorial radius [m]
      const real_t c,          // polar radius [m]
      const real_t vterm_ice,  // terminal velocity of ice particle [m/s]
      const real_t D,          // vapor diffusivity
      const real_t eta,        // air viscosity
      const real_t rho         // moist air density [kg/m^3]
      ) {

        const real_t r_max = a > c ? a : c; // particle's max radius
        const real_t X = std::cbrt(ventil::Sc(eta, rho, D)) * std::sqrt(ventil::Re(vterm_ice, r_max, rho, eta));

        const real_t b1 = (X <= real_t(1)) ? real_t(1.0)  : real_t(0.86);
        const real_t b2 = (X <= real_t(1)) ? real_t(0.14) : real_t(0.28);
        const real_t Xg = (X <= real_t(1)) ? X * X        : X; // X^gamma

        return (b1 + b2 * Xg * std::sqrt(c / capacitance(a,c)))
                / (b1 + b2 * Xg * std::sqrt(a / capacitance(a,c)));


    }
  }
}
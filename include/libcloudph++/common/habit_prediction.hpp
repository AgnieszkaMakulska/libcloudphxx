#pragma once

#include "units.hpp"
#include "const_cp.hpp"
#include "moist_air.hpp"
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
        const quantity<diffusivity, real_t> D,            // vapor diffusivity
        const quantity<thermal_conductivity, real_t> K,   // thermal conductivity
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

        return real_t(4)
          #if !defined(__NVCC__)
          * pi<real_t>()
          #else
          * CUDART_PI
          #endif
          * C * f_vnt
          * (RH_i - real_t(1))
          / (
            RH_i / D / rho_v                                   // F_d
            +
            l_s / K / T * (l_s / R_v<real_t>() / T - real_t(1)) // F_k
          );
      }

      // Capacitance of the spheroidal ice crystal
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> capacitance(
        const quantity<si::length, real_t> a,  // equatorial radius [m]
        const quantity<si::length, real_t> c  // polar radius [m]
      )
      {
        const real_t phi =  c / a;
        if (phi >= real_t(1.001)) // prolate spheroid
        {
          const real_t ecent = std::sqrt(real_t(1.0) - real_t(1.0) / (phi * phi));
          return ecent * c / std::log((real_t(1) + ecent) * phi);
        }
        if (phi <= real_t(0.999)) // oblate spheroid
        {
          const real_t ecent = std::sqrt(real_t(1) - phi * phi);
          return ecent * a / std::asin(ecent);
        }
        return (real_t(2) * a + c) / real_t(3); // approximation around phi = 1
      }


      // Inherent growth ratio Gamma(T) of Chen and Lamb (1994)
      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t growth_ratio(
        const quantity<si::temperature, real_t> T,  // temperature [K]
        const quantity<si::length, real_t> a,  // equatorial radius [m]
        const quantity<si::length, real_t> c   // polar radius [m]
      )
      {
        // quasi-spherical small crystals: Gamma = 1 if 2*max(a,c) < 10 um
        const quantity<si::length, real_t> Dmax = real_t(2) * (a > c ? a : c);
        if (Dmax < real_t(1e-5) * si::meters) return real_t(1);

        const real_t T_deg = (T - real_t(273.15) * si::kelvin) / si::kelvin;

        // nearest table entry
        int n = int(real_t(1.375) - real_t(4) * T_deg + real_t(0.5));
        if (n < 1)   n = 1;
        if (n > 121) n = 121;

        return detail::inherent_growth_ratio_table<real_t>(n - 1);
      }


      // Particle-averaged ventilation coefficient
      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t ice_ventilation_coeff(
      const quantity<si::length, real_t> a,      // equatorial radius [m]
      const quantity<si::length, real_t> c,      // polar radius [m]
      const quantity<si::velocity, real_t> vterm_ice,  // terminal velocity of ice particle [m/s]
      const quantity<diffusivity, real_t> D,          // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta,        // air viscosity
      const quantity<si::mass_density, real_t> rho         // moist air density [kg/m^3]
      ) {
        const quantity<si::length, real_t> r_max = a > c ? a : c; // particle's max radius
        const quantity<si::dimensionless, real_t> X = std::cbrt(ventil::Sc(eta, rho, D)) * std::sqrt(ventil::Re(vterm_ice, r_max, rho, eta));

        return (X <= real_t(1))
          ? real_t(1.0)  + real_t(0.14) * X * X
          : real_t(0.86) + real_t(0.28) * X;
      }


      // Primary growth habit’s ventilation coefficient
      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t ice_ventilation_coeff_ratio(
      const quantity<si::length, real_t> a,      // equatorial radius [m]
      const quantity<si::length, real_t> c,      // polar radius [m]
      const quantity<si::velocity, real_t> vterm_ice,  // terminal velocity of ice particle [m/s]
      const quantity<diffusivity, real_t> D,          // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta,        // air viscosity
      const quantity<si::mass_density, real_t> rho         // moist air density [kg/m^3]
      ) {
        const quantity<si::length, real_t> r_max = a > c ? a : c; // particle's max radius
        const quantity<si::dimensionless, real_t> X = std::cbrt(ventil::Sc(eta, rho, D)) * std::sqrt(ventil::Re(vterm_ice, r_max, rho, eta));

        const real_t b1 = (X <= real_t(1)) ? real_t(1.0)  : real_t(0.86);
        const real_t b2 = (X <= real_t(1)) ? real_t(0.14) : real_t(0.28);
        const quantity<si::dimensionless, real_t> Xg = (X <= real_t(1)) ? X * X : X; // X^gamma

        return (b1 + b2 * Xg * std::sqrt(c / capacitance(a,c)))
                / (b1 + b2 * Xg * std::sqrt(a / capacitance(a,c)));
      }

      // Deposition density
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::mass_density, real_t> rho_dep(
      const quantity<si::length, real_t> a,  // equatorial radius [m]
      const quantity<si::length, real_t> c,   // polar radius [m])
      const quantity<diffusivity, real_t> D,            // D
      const quantity<thermal_conductivity, real_t> K,   // K
      const quantity<si::mass_density, real_t> rho_v,   // ambient water vapour density
      const quantity<si::temperature, real_t> T,        // ambient temperature
      const quantity<si::dimensionless, real_t> RH_i   // p_v/p_vsi
      ) {
        if (growth_ratio(T,a,c) < real_t(1) && a < real_t(1e-4) * si::meters) return moist_air::rho_i<real_t>();

        using moist_air::R_v;
        using const_cp::p_vs, const_cp::p_vsi;
        quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t>
          l_s = const_cp::l_s<real_t>(T);

        const quantity<si::mass_density, real_t> d_rho_i =  (std::min(RH_i, p_vs(T) / p_vsi(T)) - real_t(1))
                                / D
                                / (
                                  RH_i / D / rho_v                                   // F_d
                                  +
                                  l_s / K / T * (l_s / R_v<real_t>() / T - real_t(1)) // F_k
                                );

        return  moist_air::rho_i<real_t>() * std::exp(
                      - real_t(3)
                      * std::max( (d_rho_i - real_t(5e-5) * si::kilograms / si::cubic_meters) / si::kilograms * si::cubic_meters ,real_t(0))
                      / growth_ratio(T,a,c)
                      / real_t(1e-3)
                      );
      }


      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> rho_i_new(
      const quantity<si::length, real_t> a, // equatorial radius [m]
      const quantity<si::length, real_t> c, // polar radius [m]
      const quantity<si::mass_density, real_t> rho_i, // apparent density
      const quantity<si::mass, real_t> dmi, // change of mass of an ice crystal
      const quantity<si::mass_density, real_t> rho_subldep // density of deposition / sublimation
      ) {
        quantity<si::volume, real_t> V_old = real_t(4./3.) *
                                  #if !defined(__NVCC__)
                                  * pi<real_t>()
                                  #else
                                  * CUDART_PI
                                  #endif
                                  * a * a * c;
        quantity<si::mass, real_t> mi_old = V_old * rho_i;
        quantity<si::volume, real_t> dV = dmi / rho_subldep;
        return (mi_old + dmi) / (V_old + dV);
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> ice_a_new(
      const quantity<si::length, real_t> a, // equatorial radius [m]
      const quantity<si::length, real_t> c, // polar radius [m]
      const quantity<si::mass, real_t> dmi, // change of mass of an ice crystal
      const quantity<si::mass_density, real_t> rho_subldep, // density of deposition / sublimation
      const quantity<si::temperature, real_t> T, // ambient temperature
      const quantity<si::velocity, real_t> vterm_ice,  // terminal velocity of ice particle [m/s]
      const quantity<diffusivity, real_t> D,          // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta,    // air viscosity
      const quantity<si::mass_density, real_t> rho         // moist air density [kg/m^3]
        )
      {
        quantity<si::volume, real_t> V_old = real_t(4./3.) *
                                          #if !defined(__NVCC__)
                                          * pi<real_t>()
                                          #else
                                          * CUDART_PI
                                          #endif
                                          * a * a * c;
        quantity<si::volume, real_t> V_new = V_old + dmi / rho_subldep;
        quantity<si::dimensionless, real_t> d_logV = std::log(V_new / V_old);
        real_t gamma_star;
        if (dmi > real_t(0) * si::kilograms) // deposition
        {
          gamma_star = growth_ratio(T, a, c) *  ice_ventilation_coeff_ratio(a, c, vterm_ice, D, eta, rho);
        }
        else // for sublimation gamma(T) = 1
        {
          gamma_star = ice_ventilation_coeff_ratio(a, c, vterm_ice, D, eta, rho);
        }
        return a * std::exp(d_logV / (real_t(2) + gamma_star));
      }


      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> ice_c_new(
      const quantity<si::length, real_t> a, // equatorial radius [m]
      const quantity<si::length, real_t> c, // polar radius [m]
      const quantity<si::mass, real_t> dmi, // change of mass of an ice crystal
      const quantity<si::mass_density, real_t> rho_subldep, // density of deposition / sublimation
      const quantity<si::temperature, real_t> T, // ambient temperature
      const quantity<si::velocity, real_t> vterm_ice,  // terminal velocity of ice particle [m/s]
      const quantity<diffusivity, real_t> D,          // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta,    // air viscosity
      const quantity<si::mass_density, real_t> rho         // moist air density [kg/m^3]
        )
            {
              quantity<si::volume, real_t> V_old = real_t(4./3.) *
                                                #if !defined(__NVCC__)
                                                * pi<real_t>()
                                                #else
                                                * CUDART_PI
                                                #endif
                                                * a * a * c;
              quantity<si::volume, real_t> V_new = V_old + dmi / rho_subldep;
              quantity<si::dimensionless, real_t> d_logV = std::log(V_new / V_old);
              real_t gamma_star;
              if (dmi > real_t(0) * si::kilograms) // deposition
              {
                gamma_star = growth_ratio(T, a, c) *  ice_ventilation_coeff_ratio(a, c, vterm_ice, D, eta, rho);
              }
              else // for sublimation gamma(T) = 1
              {
                gamma_star = ice_ventilation_coeff_ratio(a, c, vterm_ice, D, eta, rho);
              }
              return c * std::exp(gamma_star * d_logV / (real_t(2) + gamma_star));
      }


    }
  }
}
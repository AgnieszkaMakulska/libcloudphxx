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

      // Capacitance of the spheroidal ice crystal
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> capacitance(
        const quantity<si::length, real_t> a, // equatorial radius
        const quantity<si::length, real_t> c  // polar radius
      )
      {
        assert(a > real_t(0) * si::meters);
        assert(c > real_t(0) * si::meters);
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
      quantity<si::dimensionless, real_t> growth_ratio(
        const quantity<si::temperature, real_t> T, // temperature
        const quantity<si::length, real_t> a,      // equatorial radius
        const quantity<si::length, real_t> c       // polar radius
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
      quantity<si::dimensionless, real_t> ice_ventilation_coeff(
      const quantity<si::length, real_t> a,              // equatorial radius
      const quantity<si::length, real_t> c,              // polar radius
      const quantity<si::velocity, real_t> vterm_ice,    // terminal velocity of ice particle
      const quantity<diffusivity, real_t> D,             // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta, // air viscosity
      const quantity<si::mass_density, real_t> rho       // moist air density
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
      quantity<si::dimensionless, real_t> ice_ventilation_coeff_ratio(
      const quantity<si::length, real_t> a,              // equatorial radius
      const quantity<si::length, real_t> c,              // polar radius
      const quantity<si::velocity, real_t> vterm_ice,    // terminal velocity of ice particle
      const quantity<diffusivity, real_t> D,             // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta, // air viscosity
      const quantity<si::mass_density, real_t> rho       // moist air density
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
      const quantity<si::length, real_t> a,           // equatorial radius
      const quantity<si::length, real_t> c,           // polar radius
      const quantity<diffusivity, real_t> D,          // vapor diffusivity
      const quantity<thermal_conductivity, real_t> K, // thermal conductivity
      const quantity<si::mass_density, real_t> rho_v, // ambient water vapour density
      const quantity<si::temperature, real_t> T,      // ambient temperature
      const quantity<si::dimensionless, real_t> RH_i  // p_v/p_vsi
      ) {
        if (growth_ratio(T,a,c) < real_t(1) && a < real_t(1e-4) * si::meters) return moist_air::rho_i_true<real_t>();

        using moist_air::R_v;
        using const_cp::p_vs, const_cp::p_vsi;
        quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t>
          l_s = const_cp::l_s<real_t>(T);

        const quantity<si::dimensionless, real_t> RH_lim = RH_i < p_vs(T) / p_vsi(T) ? RH_i : p_vs(T) / p_vsi(T);
        const real_t d_rho_i =  (RH_lim - real_t(1))
                                / D
                                / (
                                  RH_i / D / rho_v                                   // F_d
                                  +
                                  l_s / K / T * (l_s / R_v<real_t>() / T - real_t(1)) // F_k
                                ) / si::kilograms * si::cubic_meters;
        const real_t d_rho_max =
                      (d_rho_i - real_t(5e-5)) > real_t(0) ?
                       d_rho_i - real_t(5e-5) : real_t(0);

        return  moist_air::rho_i_true<real_t>() * std::exp(
                      - real_t(3) * d_rho_max
                      / growth_ratio(T,a,c)
                      / real_t(1e-3)
                      );
      }


      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::mass_density, real_t> ice_rho_update(
      const quantity<si::length, real_t> a,                // equatorial radius
      const quantity<si::length, real_t> c,                // polar radius
      const quantity<si::mass_density, real_t> rho_i,      // apparent density
      const quantity<si::mass, real_t> dmi,                // change of mass of an ice crystal
      const quantity<si::mass_density, real_t> rho_subldep // density of deposition / sublimation
      ) {
        assert(a > real_t(0) * si::meters);
        assert(c > real_t(0) * si::meters);
        quantity<si::volume, real_t> V_old = real_t(4./3.)
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
      quantity<si::length, real_t> ice_a_update(
      const quantity<si::length, real_t> a,                 // equatorial radius
      const quantity<si::length, real_t> c,                 // polar radius
      const quantity<si::mass, real_t> dmi,                 // change of mass of an ice crystal
      const quantity<si::mass_density, real_t> rho_subldep, // density of deposition / sublimation
      const quantity<si::temperature, real_t> T,            // ambient temperature
      const quantity<si::velocity, real_t> vterm_ice,       // terminal velocity of ice particle
      const quantity<diffusivity, real_t> D,                // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta,    // air viscosity
      const quantity<si::mass_density, real_t> rho          // moist air density
        )
      {
        assert(a > real_t(0) * si::meters);
        assert(c > real_t(0) * si::meters);
        quantity<si::volume, real_t> V_old = real_t(4./3.)
                                          #if !defined(__NVCC__)
                                          * pi<real_t>()
                                          #else
                                          * CUDART_PI
                                          #endif
                                          * a * a * c;
        quantity<si::volume, real_t> V_new = V_old + dmi / rho_subldep;
        quantity<si::dimensionless, real_t> d_logV = std::log(V_new / V_old);
        real_t gamma_star;
        if (dmi >= real_t(0) * si::kilograms) // deposition
        {
          gamma_star = real_t(c/a) < real_t(40) ?   // limiting the growth ratio of very thin crystals
            real_t(growth_ratio(T, a, c) *  ice_ventilation_coeff_ratio(a, c, vterm_ice, D, eta, rho))
            : real_t(1);
        }
        else // for sublimation gamma(T) = 1
        {
          gamma_star = ice_ventilation_coeff_ratio(a, c, vterm_ice, D, eta, rho);
        }
        return a * std::exp(d_logV / (real_t(2) + gamma_star));
      }


      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> ice_c_update(
      const quantity<si::length, real_t> a,                 // equatorial radius
      const quantity<si::length, real_t> c,                 // polar radius
      const quantity<si::mass, real_t> dmi,                 // change of mass of an ice crystal
      const quantity<si::mass_density, real_t> rho_subldep, // density of deposition / sublimation
      const quantity<si::temperature, real_t> T,            // ambient temperature
      const quantity<si::velocity, real_t> vterm_ice,       // terminal velocity of ice particle
      const quantity<diffusivity, real_t> D,                // vapor diffusivity
      const quantity<si::dynamic_viscosity, real_t> eta,    // air viscosity
      const quantity<si::mass_density, real_t> rho          // moist air density
        )
            {
              assert(a > real_t(0) * si::meters);
              assert(c > real_t(0) * si::meters);
              quantity<si::volume, real_t> V_old = real_t(4./3.)
                                                #if !defined(__NVCC__)
                                                * pi<real_t>()
                                                #else
                                                * CUDART_PI
                                                #endif
                                                * a * a * c;
              quantity<si::volume, real_t> V_new = V_old + dmi / rho_subldep;
              quantity<si::dimensionless, real_t> d_logV = std::log(V_new / V_old);
              real_t gamma_star;
              if (dmi >= real_t(0) * si::kilograms) // deposition
              {
                gamma_star = real_t(c/a) < real_t(40) ?   // limiting the growth ratio of very thin crystals
                  real_t(growth_ratio(T, a, c) *  ice_ventilation_coeff_ratio(a, c, vterm_ice, D, eta, rho))
                  : real_t(1);
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
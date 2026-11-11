#pragma once

#include "units.hpp"
#include "macros.hpp"
#include "const_cp.hpp"

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace common
  {
    namespace maxwell_mason
    {
      // Condensational growth rate of a liquid droplet
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::area, si::time>::type, real_t> rdrdt(
        const quantity<diffusivity, real_t> D,            // D 
        const quantity<thermal_conductivity, real_t> K,   // K
        const quantity<si::mass_density, real_t> rho_v,   // ambient water vapour density
        const quantity<si::temperature, real_t> T,        // ambient temperature
        const quantity<si::pressure, real_t> p,           // ambient pressure
        const quantity<si::dimensionless, real_t> RH,     // p_v/p_vs = relative humidity
        const quantity<si::dimensionless, real_t> a_w,    // water activity
        const quantity<si::dimensionless, real_t> klvntrm // the Kelvin term
      )
      {
        using moist_air::rho_w;
        using moist_air::R_v;

        quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t> 
          l_v = const_cp::l_v<real_t>(T);

        return (real_t(1) - a_w * klvntrm / RH)
          / rho_w<real_t>() 
          / ( 
            real_t(1) 
              / D
              / rho_v 
            +
            l_v 
              / K
              / RH  
              / T 
              * (l_v / R_v<real_t>() / T - real_t(1))
          )   
        ;   
      }

      // Depositional growth rate of ice: dm/dt following Shima et al. 2020
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dmi_dt(
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
    };
  };
};

// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <cmath>

// #include <thrust/iterator/transform_iterator.h>
#include <libcloudph++/common/habit_prediction.hpp>
#include <libcloudph++/common/maxwell-mason.hpp>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // Mass-growth functor. dmi_dt(mi) is the rate evaluated for a trial mass mi:
      // the trial geometry (a, c) is obtained. operator() is the residual of the backward
      // Euler scheme (to be used with toms748); the explicit scheme uses dmi_dt()
      // evaluated at the old mass.

      template <typename real_t>
      struct advance_ice_mass_minfun
      {
        const quantity<si::time,              real_t> dt;
        const quantity<si::length,            real_t> a_old;
        const quantity<si::length,            real_t> c_old;
        const quantity<si::mass_density,      real_t> rho_i_old;
        const quantity<si::mass_density,      real_t> rhod;
        const quantity<si::dimensionless,     real_t> rv;
        const quantity<si::temperature,       real_t> T;
        const quantity<si::dynamic_viscosity, real_t> eta;
        const quantity<si::velocity,          real_t> vt;
        const quantity<si::pressure,          real_t> p;
        const quantity<si::dimensionless,     real_t> RH_i;
        const quantity<si::dimensionless,     real_t> RH_max;

        quantity<si::mass, real_t> mi_old;              // old ice crystal mass
        quantity<si::mass_density, real_t> rho_subldep; // density of deposition or sublimation

        // ctor
        BOOST_GPU_ENABLED
        advance_ice_mass_minfun(
          const real_t &dt,
          const real_t &a,
          const real_t &c,
          const real_t &rho_i,
          const thrust::tuple<thrust::tuple<real_t, real_t, real_t, real_t, real_t>, real_t, real_t> &tpl,
          const real_t &RH_max
        ) :
          dt(dt * si::seconds),
          a_old(a * si::meters),
          c_old(c * si::meters),
          rho_i_old(rho_i * si::kilograms / si::cubic_meters),
          rhod(    thrust::get<0>(thrust::get<0>(tpl)) * si::kilograms / si::cubic_metres),
          rv(      thrust::get<1>(thrust::get<0>(tpl))),
          T(       thrust::get<2>(thrust::get<0>(tpl)) * si::kelvins),
          eta(     thrust::get<3>(thrust::get<0>(tpl)) * si::pascals * si::seconds),
          vt(      thrust::get<4>(thrust::get<0>(tpl)) * si::metres_per_second),
          p(       thrust::get<1>(tpl) * si::pascals),
          RH_i(    thrust::get<2>(tpl)),
          RH_max(RH_max)
        {
          using common::habit_prediction::rho_dep;
          using common::moist_air::D_0;
          using common::moist_air::K_0;

          #if !defined(__NVCC__)
                    const real_t pi_ = pi<real_t>();
          #else
                    const real_t pi_ = CUDART_PI;
          #endif

          mi_old = rho_i_old * real_t(4./3.) * pi_ * a_old * a_old * c_old;

          // the sign of dm/dt equals the sign of (RH_i - 1)
          if (RH_i >= real_t(1))
            rho_subldep = rho_dep<real_t>(a_old, c_old, D_0<real_t>(), K_0<real_t>(), rhod * rv, T, RH_i);
          else
            rho_subldep = rho_i_old; // sublimation does not change the apparent density
        }

        // dm/dt for a trial mass
        BOOST_GPU_ENABLED
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t>
        dmi_dt(const quantity<si::mass, real_t> &mi_in) const
        {
          using namespace common::habit_prediction;
          using common::moist_air::D_0;
          using common::moist_air::K_0;

          // prohibit complete sublimation
          const quantity<si::mass, real_t> mi    = mi_in > real_t(3.84e-24) * si::kilograms ?
                                                   mi_in : real_t(3.84e-24) * si::kilograms;
          const quantity<si::mass, real_t> dmi   = mi - mi_old;

          // no ventilation or kinetic correction in D and K, ventilation is defined by f_vnt
          const quantity<common::diffusivity, real_t>          D = D_0<real_t>();
          const quantity<common::thermal_conductivity, real_t> K = K_0<real_t>();

          // geometry at the trial mass
          // (if mi_in = m_old, then a, c = a_old, c c_old)
          const quantity<si::length, real_t> a =
            ice_a_update<real_t>(a_old, c_old, dmi, rho_subldep, T, vt, D, eta, rhod);
          const quantity<si::length, real_t> c =
            ice_c_update<real_t>(a_old, c_old, dmi, rho_subldep, T, vt, D, eta, rhod);

          const quantity<si::length, real_t>        C     = capacitance<real_t>(a, c);
          const quantity<si::dimensionless, real_t> f_vnt = ice_ventilation_coeff<real_t>(a, c, vt, D, eta, rhod);

          return common::maxwell_mason::dmi_dt<real_t>(
            D,          // vapor diffusivity
            K,          // thermal conductivity
            rhod * rv,  // ambient water vapour density
            T,          // ambient temperature
            RH_i,       // p_v/p_vsi
            C,          // capacitance of the crystal
            f_vnt       // ventilation coeff.
          );
        }

        // backward Euler scheme: minimize dm/dt(mi) * dt + mi_old - mi
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &mi_unitless) const
        {
          const quantity<si::mass, real_t> mi = mi_unitless * si::kilograms;
          return (mi_old + dt * dmi_dt(mi) - mi) / si::kilograms;
        }
      };


      template <typename real_t, bool apply = true>
      struct advance_ice_attributes
      {
        const real_t dt, RH_max;

        BOOST_GPU_ENABLED
        advance_ice_attributes(const real_t &dt, const real_t &RH_max)
          : dt(dt), RH_max(RH_max) {}

        BOOST_GPU_ENABLED
        thrust::tuple<real_t, real_t, real_t> operator()(
            const thrust::tuple<real_t, real_t, real_t> &ice_attributes_old,
            const thrust::tuple<
                thrust::tuple<real_t, real_t, real_t, real_t, real_t>,
                real_t, real_t> &tpl
        ) const
        {
          using namespace common::habit_prediction;
          using common::moist_air::D_0;
          using common::moist_air::K_0;
          using common::moist_air::rho_i;

          #if !defined(__NVCC__)
                    const real_t pi_ = pi<real_t>();
          #else
                    const real_t pi_ = CUDART_PI;
          #endif

          const real_t a_old     = thrust::get<0>(ice_attributes_old);
          const real_t c_old     = thrust::get<1>(ice_attributes_old);
          const real_t rho_i_old = thrust::get<2>(ice_attributes_old);

          // skip liquid droplets
          if (a_old <= real_t(0) || c_old <= real_t(0))
          {
            if constexpr (apply)
              return ice_attributes_old;
            else
              return thrust::make_tuple(real_t(0), real_t(0), real_t(0));
          }

          const thrust::tuple<real_t, real_t, real_t, real_t, real_t> &st = thrust::get<0>(tpl);
          const quantity<si::mass_density, real_t>      rhod = thrust::get<0>(st) * si::kilograms / si::cubic_metres;
          const quantity<si::dimensionless, real_t>     rv   = thrust::get<1>(st);
          const quantity<si::temperature, real_t>       T    = thrust::get<2>(st) * si::kelvins;
          const quantity<si::dynamic_viscosity, real_t> eta  = thrust::get<3>(st) * si::pascals * si::seconds;
          const quantity<si::velocity, real_t>          vt   = thrust::get<4>(st) * si::metres_per_second;


          const quantity<common::diffusivity, real_t> D = D_0<real_t>();

          const advance_ice_mass_minfun<real_t> f(dt, a_old, c_old, rho_i_old, tpl, RH_max);

          const real_t mi_old = f.mi_old / si::kilograms;

          // Explicit (forward) Euler using the rate at the old state.
          const real_t dmi_dt_val = f.dmi_dt(f.mi_old) * si::seconds / si::kilograms; // kg/s as a number
          real_t dmi = dt * dmi_dt_val;

          // TODO: implicit integration

          const real_t m_min = real_t(3.84e-24);
          if (dmi < m_min - mi_old) dmi = m_min - mi_old;
          real_t mi_new = mi_old + dmi;

          // multiplying by units
          const quantity<si::length, real_t> a_q = a_old * si::meters;
          const quantity<si::length, real_t> c_q = c_old * si::meters;
          const quantity<si::mass, real_t>   dmi_q = dmi * si::kilograms;
          const quantity<si::mass_density, real_t> rho_i_old_q = rho_i_old * si::kilograms / si::cubic_meters;

          real_t a_new = ice_a_update<real_t>(a_q, c_q, dmi_q, f.rho_subldep, T, vt, D, eta, rhod) / si::meters;
          real_t c_new = ice_c_update<real_t>(a_q, c_q, dmi_q, f.rho_subldep, T, vt, D, eta, rhod) / si::meters;
          real_t rho_i_new = ice_rho_update<real_t>(a_q, c_q, rho_i_old_q, dmi_q, f.rho_subldep)
                             / (si::kilograms / si::cubic_meters);

          // if the minor axis < 1 um, assume it's a sphere with the true ice density
          if ((a_new < c_new ? a_new : c_new) < real_t(1e-6))
          {
            rho_i_new = rho_i<real_t>() / (si::kilograms / si::cubic_meters);
            a_new = c_new = std::cbrt(mi_new / (rho_i_new * real_t(4./3.) * pi_));
          }

          if constexpr (apply)
            return thrust::make_tuple(a_new, c_new, rho_i_new);
          else
            return thrust::make_tuple(a_new - a_old, c_new - c_old, rho_i_new - rho_i_old);
        }
      };

    };
  };
};

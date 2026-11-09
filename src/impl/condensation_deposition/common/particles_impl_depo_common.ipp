// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// #include <thrust/iterator/transform_iterator.h>
#include <libcloudph++/common/habit_prediction.hpp>
#include <libcloudph++/common/maxwell-mason.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>
#include <libcloudph++/common/kelvin_term.hpp>
#include <libcloudph++/common/transition_regime.hpp>
#include <libcloudph++/common/ventil.hpp>
#include <libcloudph++/common/mean_free_path.hpp>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct advance_ice_mass_minfun
      {
        const quantity<si::length,            real_t> a;
        const quantity<si::length,            real_t> c;
        const quantity<si::mass_density,      real_t> rho_i;
        const quantity<si::time,              real_t> dt;
        const quantity<si::mass_density,      real_t> rhod;
        const quantity<si::dimensionless,     real_t> rv;
        const quantity<si::temperature,       real_t> T;
        const quantity<si::pressure,          real_t> p;
        const quantity<si::dimensionless,     real_t> RH_i;
        const quantity<si::dynamic_viscosity, real_t> eta;
        const quantity<si::velocity,          real_t> vt;
        const quantity<si::dimensionless,     real_t> RH_max;
        const quantity<si::length,            real_t> lambda_D;
        const quantity<si::length,            real_t> lambda_K;

        // ctor
        BOOST_GPU_ENABLED
        advance_ice_mass_minfun(
          const real_t &dt,
          const real_t &a,
          const real_t &c,
          const real_t &rho_i,
          const thrust::tuple<thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t>, real_t, real_t> &tpl,
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
          RH_i(      thrust::get<2>(tpl)),
          lambda_D(thrust::get<5>(thrust::get<0>(tpl)) * si::metres),
          lambda_K(thrust::get<6>(thrust::get<0>(tpl)) * si::metres),
          RH_max(RH_max)
        {}

        BOOST_GPU_ENABLED
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dmi_dt(const quantity<si::mass, real_t> &mi) const
        {
          using namespace common::habit_prediction;
          using namespace common::maxwell_mason;
          using namespace common::kappa_koehler;
          using namespace common::kelvin;
          using common::moist_air::D_0;
          using common::moist_air::K_0;
          using common::moist_air::c_pd;
          using common::transition_regime::beta;
          using common::ventil::Sh;
          using common::ventil::Nu;
#if !defined(__NVCC__)
          using std::sqrt;
#endif

          const quantity<si::dimensionless, real_t>
            Re = common::ventil::Re(vt, rw, rhod, eta),
            Sc = common::ventil::Sc(eta, rhod, D_0<real_t>()), // TODO? cache
            Pr = common::ventil::Pr(eta, c_pd<real_t>(), K_0<real_t>()); // TODO? cache

          const quantity<common::diffusivity, real_t>
            D = D_0<real_t>() * beta(lambda_D / rw) * (Sh(Sc, Re) / 2);

          const quantity<common::thermal_conductivity, real_t>
            K = K_0<real_t>() * beta(lambda_K / rw) * (Nu(Pr, Re) / 2);

          return dmi_dt(

            )
        }

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &mi_unitless) const
        {
          const quantity<si::mass, real_t> mi = mi_unitless * si::square_metres;
          return (mi_old + dt * dmi_dt(mi) - mi) / si::square_metres;
        }
      };




      template <typename real_t, bool apply = true>
      struct advance_ice_attributes
      {
        const real_t dt, RH_max;
        const quantity<si::temperature,              real_t> T;
        const quantity<si::dimensionless,            real_t> RH_i;
        const quantity<si::mass_density,             real_t> rho;
        const quantity<si::dynamic_viscosity,        real_t> eta;
        const quantity<common::diffusivity,          real_t> D;
        const quantity<common::thermal_conductivity, real_t> K;
        const quantity<si::velocity,                 real_t> vterm_ice;

        advance_ice_attributes(const real_t &dt, const real_t &RH_max)
          : dt(dt), RH_max(RH_max) {}

        BOOST_GPU_ENABLED
        thrust::tuple<real_t, real_t> operator()(
            const thrust::tuple<real_t, real_t> &ice_attributes_old,
            const thrust::tuple<
                thrust::tuple<real_t, real_t, real_t, real_t, real_t, real_t, real_t>,
                real_t, real_t> &tpl
        ) const
        {
          const real_t a_old = thrust::get<0>(ice_attributes_old);
          const real_t c_old = thrust::get<1>(ice_attributes_old);
          const real_t rho_i_old = thrust::get<2>(ice_attributes_old);
          const real_t mi_old = rho_i_old * real_t(4./3.)
                                #if !defined(__NVCC__)
                                * pi<real_t>()
                                #else
                                * CUDART_PI
                                #endif
                                * a_old * a_old * c_old;

          // Skip liquid droplets
          if (a_old <= 0 || c_old <= 0)       
          {
            if constexpr (apply)
              return ice_attributes_old;
            else
              return thrust::make_tuple(real_t(0), real_t(0), real_t(0));
          }

          // CALCULATE THE NEW MASS
          // forward Euler for simplicity, TODO: implicit Euler / toms748 as in condensation
          // const real_t dm_dt =

          // limit dmi to prohibit complete sublimation
          dmi = std::max(dmi, real_t(4e-24) - mi_old);

          real_t rho_subldep, a_new, c_new, rho_i_new;

          if (dmi >= real_t(0))  // deposition
          {
            rho_subldep = common::habit_prediction::rho_dep(a_old, c_old, D, K, rhod * rv, T, RH_i);
          }
          else  // sublimation does not change the density
          {
            rho_subldep = rho_i_old;
          }

          a_new = common::habit_prediction::ice_a_update<real_t>(a_old, c_old, dmi, rho_subldep, T, vterm_ice, D, eta, rho);
          c_new = common::habit_prediction::ice_c_update<real_t>(a_old, c_old, dmi, rho_subldep, T, vterm_ice, D, eta, rho);
          rho_i_new = common::habit_prediction::ice_rho_update<real_t>(a_old, c_old, rho_i_old, dmi, rho_subldep);

          // if the minor axis < 1 um, assume it's spherical with true ice density
          if ((a_new < c_new ? a_new : c_new) < real_t(1e-6) * si::metres)
          {
            rho_i_new = common::moist_air::rho_i<real_t>();
            a_new = std::cbrt(mi_new / (rho_i_new * real_t(4./3.)
                            #if !defined(__NVCC__)
                            * pi<real_t>()
                            #else
                            * CUDART_PI
                            #endif
                            )) * si::metres;
            c_new = std::cbrt(mi_new / (rho_i_new * real_t(4./3.)
                            #if !defined(__NVCC__)
                            * pi<real_t>()
                            #else
                            * CUDART_PI
                            #endif
                            )) * si::metres;
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

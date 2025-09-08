// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
   * @brief GPU-accelerated routines for super-droplet condensation and 3rd wet moment updates.
   *
   * Handles:
   * - Sorting particles by grid cell
   * - Computing pre- and post-condensation 3rd wet moments
   * - Updating particle radii via backward-Euler integration
   * - Accumulating net 3rd moment changes
   * - Updating Eulerian thermodynamic fields (rv, th)
   *
   * All operations use Thrust vectors and zip/permutation iterators to
   * couple Lagrangian particles with Eulerian grid fields.
   */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::cond(
      const real_t &dt,
      const real_t &RH_max,
      const bool turb_cond
    ) {   

      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &lambda_D(tmp_device_real_cell1); // real_cell used in cond.ipp
      thrust_device::vector<real_t> &lambda_K(tmp_device_real_cell2); // real_cell used in cond.ipp

      // calculating liquid water content before condensation
      hskpng_sort(); // sorting particles by cell index
      thrust_device::vector<real_t> &drv(tmp_device_real_cell);

      // calculating the 3rd wet moment before condensation from rw2
      moms_all(); //prepare particles for moment calculation
      moms_calc(rw2.begin(), real_t(3./2.)); // calculate the 3rd moment
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) before condensation");

      // fill drv with 0s if not all cells will be updated in the following transform
      if(count_n!=n_cell)  thrust::fill(drv.begin(), drv.end(), real_t(0.));

      // permute-copying the result to -dm_3
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - the 3rd moment before condensation
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output - drv at indices given by count_ijk
        thrust::negate<real_t>()                                               // multiplying by -1
      );

      // tuple (rhod, rv, T, eta, rd3, kpa, vt, lambda_D, lambda_K) for each particle
      auto hlpr_zip_iter = thrust::make_zip_iterator(thrust::make_tuple(
        thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
        thrust::make_permutation_iterator(rv.begin(), ijk.begin()),
        thrust::make_permutation_iterator(T.begin(), ijk.begin()),
        thrust::make_permutation_iterator(eta.begin(), ijk.begin()),
        rd3.begin(),
        kpa.begin(),
        vt.begin(),
        thrust::make_permutation_iterator(lambda_D.begin(), ijk.begin()),
        thrust::make_permutation_iterator(lambda_K.begin(), ijk.begin())
      ));

      // calculating drop growth in a timestep using backward Euler 
      // TODO: both calls almost identical, use std::bind or sth?
      if(turb_cond)
      {
        // if turb_cond, adding the particle supersaturation contribution to the interpolated grid RH
        thrust_device::vector<real_t> &RH_plus_ssp(tmp_device_real_part2);
        thrust::transform(
          ssp.begin(), ssp.end(),
          thrust::make_permutation_iterator(RH.begin(), ijk.begin()),
          RH_plus_ssp.begin(),
          arg::_1 + arg::_2
        );

        thrust::transform(
          rw2.begin(), rw2.end(),         // input 1: current rw2 values
          thrust::make_zip_iterator(      // input 2: tuple of particle & cell info
            thrust::make_tuple(
              hlpr_zip_iter,               // particle+cell state (rhod, rv, T, eta, rd3, ...)
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              RH_plus_ssp.begin()              // total supersaturation for each particle
            )
          ), 
          rw2.begin(),                    // output: updated rw2
          detail::advance_rw2<real_t>(dt, RH_max) // functor that does the Euler step
        );
      }
      else
        thrust::transform(
          rw2.begin(), rw2.end(),         // input 1: current rw2 values
          thrust::make_zip_iterator(      // input 2: tuple of particle & cell info
            thrust::make_tuple(
              hlpr_zip_iter,               // particle+cell state (rhod, rv, T, eta, rd3, ...)
              thrust::make_permutation_iterator(p.begin(), ijk.begin()),
              thrust::make_permutation_iterator(RH.begin(), ijk.begin()) // cell RH only
            )
          ), 
          rw2.begin(),                    // output: updated rw2
          detail::advance_rw2<real_t>(dt, RH_max) // functor that does the Euler step
        );
      nancheck(rw2, "rw2 after condensation (no sub-steps");

      // calculating the 3rd wet moment after condensation from rw2
      moms_calc(rw2.begin(), real_t(3./2.));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) after condensation");

      // adding the 3rd moment after condensation to the existing drv
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - the 3rd moment after condensation
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // input - count_mom values mapped to cell indices in drv
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()                                                // adding
      );

      // update th and rv according to changes in third specific wet moment
      update_th_rv(drv);
    }
  };  
};

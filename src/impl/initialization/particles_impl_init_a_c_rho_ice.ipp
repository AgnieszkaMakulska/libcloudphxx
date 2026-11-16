// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
    namespace lgrngn
    {
        template <typename real_t, backend_t device>
        void particles_t<real_t, device>::impl::init_a_c_rho_ice()
        {
            if (opts_init.init_ice)
            {
            // initialize spherical ice crystals of radius ice_r0 and the true ice density
            thrust::fill(ice_a.begin() + n_part_old, ice_a.end(), opts_init.ice_r0);
            thrust::fill(ice_c.begin() + n_part_old, ice_c.end(), opts_init.ice_r0);
            thrust::fill(ice_rho.begin() + n_part_old, ice_rho.end(),
              real_t(common::moist_air::rho_i_true<real_t>() / si::kilograms * si::cubic_metres));
            }
            else
            {
            // filling a, c, rho_i with zeros initially
            thrust::fill(ice_a.begin() + n_part_old, ice_a.end(), real_t(0));
            thrust::fill(ice_c.begin() + n_part_old, ice_c.end(), real_t(0));
            thrust::fill(ice_rho.begin() + n_part_old, ice_rho.end(), real_t(0));
            }
        }
    };
};
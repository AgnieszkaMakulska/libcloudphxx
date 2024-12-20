#pragma once

#include "spectral.hpp"

template <typename real_t>
class settings_t {
public:
    const std::string aerosol, init;
    const real_t
            z_max = 200, // [m]
            vertical_velocity = 2.,
            dt = 0.1,
            kappa = 0.61;
    real_t thd0 =292.; // [K]
    real_t rv0 = 0.008; // [kg/kg]
    real_t rhod0 = 1.1; // [kg/m3]
    const int sd_conc = 2000;
    const int n_sd_max = 2000;
    const int n_steps = z_max / (vertical_velocity * dt);

    std::shared_ptr<bimodal<real_t>> n_ln_rd_stp =  std::make_shared<bimodal<real_t>>(
              lognormal<real_t>(0.03e-6, 1.28, 1000.0e6),
              lognormal<real_t>(0.14e-6, 1.75, 50.0e6)
      );
};
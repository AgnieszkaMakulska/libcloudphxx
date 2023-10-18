#include <ranges>

#include "libcloud_hacks.hpp"
#include "output.hpp"
#include "spectral.hpp"

#include <libcloudph++/lgrngn/factory.hpp>
#include <libcloudph++/lgrngn/terminal_velocity.hpp>
#include <libcloudph++/lgrngn/kernel.hpp>

using real_t = double;
const auto backend = libcloudphxx::lgrngn::serial;


int main(int arg_count, char** arg_values)
{
    libcloudphxx::lgrngn::opts_init_t<real_t> params;
    libcloudphxx::lgrngn::opts_t<real_t> opts;

    real_t r_zero = 30.084e-6;
    real_t n_zero = 1.25 * pow(2,23);

    real_t dv = 100; // [m3]
    params.dx = params.dy = params.dz = pow(dv, real_t(1./3));
    params.nx = params.ny = params.nz = 0;
    params.dt = 1;
    //params.n_sd_max = pow(10.,5);
    //params.sd_const_multi = n_zero/params.n_sd_max;
    params.sd_conc =  pow(10.,5);
    params.n_sd_max = params.sd_conc;

    real_t kappa=1e-10;
    std::shared_ptr<exponential<real_t>> n_ln_rd = std::make_shared<exponential<real_t>>(r_zero, n_zero);
    params.dry_distros.emplace(kappa, n_ln_rd);

    params.kernel = libcloudphxx::lgrngn::kernel_t::hall_davis_no_waals;
    params.sstp_coal=1;
    params.terminal_velocity=libcloudphxx::lgrngn::vt_t::beard76;


    opts.adve = false;
    params.sedi_switch = opts.sedi = false;
    opts.cond = false;
    opts.coal = true;

    std::unique_ptr<
            libcloudphxx::lgrngn::particles_proto_t<real_t>
    > prtcls(libcloudphxx::lgrngn::factory<real_t>(
            backend,
            params
    ));

    real_t
            thd = 300.,
            rv = 0.01,
            rhod = 1.;

    real_t simulation_time=1800.;
    auto n_steps = int(simulation_time / params.dt);

    std::string filename = "data.nc";
    auto range_i = std::views::iota(0, n_steps + 1);
    auto nc = output_init(int(params.n_sd_max), int(range_i.size()), params.dt, filename);

    prtcls->init(arrinfo(thd), arrinfo(rv), arrinfo(rhod));
    {
        for (auto i: range_i) {
            if (i != 0) {
                prtcls->step_sync(opts, arrinfo(thd), arrinfo(rv), arrinfo(rhod));
                output_step<backend>(i, *prtcls, *nc);
                prtcls->step_async(opts);
            }
            else
                output_step<backend>(i, *prtcls, *nc);
        }
    }
};
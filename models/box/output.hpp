#pragma once
#include <netcdf>
#include <valarray>
#include <boost/range/combine.hpp>
#include "libcloud_hacks.hpp"
#include <cmath>
#include <vector>

int numBins = 150; //number of bins for spectra

template <typename real_t>
auto output_init(
        int n_sd,
        int n_t,
        real_t dt,
        std::string &filename
) {
    auto nc = std::make_shared<netCDF::NcFile>(filename, netCDF::NcFile::replace);

    nc->putAtt("dt (s)", netCDF::ncFloat, dt);

    nc->addDim("step",n_t);
    nc->addDim("droplet_id",n_sd);
    nc->addDim("numBin", numBins);

    nc->addVar("wet radius squared", "double", std::vector<std::string>{"step", "droplet_id"}).putAtt("unit", "m^2");
    nc->addVar("mass density", "double", std::vector<std::string>{"step", "numBin"}).putAtt("unit", "g/ (m^3 lnr)");
    nc->addVar("moment 0", "double", std::vector<std::string>{"step", "numBin"}).putAtt("unit", "1/kg");
    nc->addVar("moment 3", "double", std::vector<std::string>{"step", "numBin"}).putAtt("unit", "m^3/kg");
    return nc;
}

template <typename real_t>
auto save_scalar(
        const int i,
        libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls,
        const netCDF::NcFile &nc,
        const std::string &name
) {
    auto value = prtcls.outbuf()[0];
    nc.getVar(name).putVar(std::vector{size_t(i)}, value);
    return value;
}


template <typename iter_t>
void save_vector(
        const int i,
        iter_t &iter,
        const netCDF::NcFile &nc,
        const std::string &name
)
{
    for (size_t j=0; const auto &item : iter)
        nc.getVar(name).putVar(std::vector{size_t(i), j++}, item);
}


template <
        libcloudphxx::lgrngn::backend_t backend,
        typename real_t
>
void output_step(
        const int i,
        libcloudphxx::lgrngn::particles_proto_t<real_t> &prtcls,
        const netCDF::NcFile &nc
) {


    //save_vector(i, impl<backend, real_t>(prtcls)->rw2, nc, "wet radius squared");
    std::vector<real_t> rw2 = prtcls.get_attr("rw2");
    save_vector(i, rw2, nc, "wet radius squared");

    std::vector<real_t> bins(numBins);
    for (int j = 0; j < numBins; j++) {
        bins[j] = 6.0 * std::pow(10, -6 + j / 50.0);
    }

    std::vector<real_t> mass_density(numBins);
    for (int j = 0; j < numBins; j++) {
        prtcls.diag_all();
        prtcls.diag_wet_mass_dens( (bins[j] + bins[j+1])/2. ,0.62 );
        auto value = prtcls.outbuf()[0];
        mass_density[j] = value * 1000; // *1000 to get grams
    }
    save_vector(i, mass_density, nc, "mass density");

    std::vector<real_t> mom_0(numBins);
    for (int j = 0; j < numBins; j++) {
        prtcls.diag_wet_rng(bins[j], bins[j+1]);
        prtcls.diag_wet_mom(0);
        auto value = prtcls.outbuf()[0];
        mom_0[j] = value;
    }
    save_vector(i, mom_0, nc, "moment 0");

    std::vector<real_t> mom_3(numBins);
    for (int j = 0; j < numBins; j++) {
        prtcls.diag_wet_rng(bins[j], bins[j+1]);
        prtcls.diag_wet_mom(3);
        auto value = prtcls.outbuf()[0];
        mom_3[j] = value;
    }
    save_vector(i, mom_3, nc, "moment 3");


};
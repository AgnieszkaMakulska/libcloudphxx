#pragma once
#include <libcloudph++/common/unary_function.hpp>

template <typename real_t>
struct exponential : libcloudphxx::common::unary_function<real_t> {
    const real_t r_zero, n_zero;

    explicit exponential(
            const real_t r_zero, const real_t n_zero
    ) : r_zero(r_zero), n_zero(n_zero)
    {}

    real_t funval(const real_t lnr) const {
        real_t r = exp(lnr);
        return n_zero * 3.*pow(r,3) / pow(r_zero,3) * exp(- pow((r/r_zero),3));
    }
};

template <typename real_t>
struct bimodal : libcloudphxx::common::unary_function<real_t> {
    const exponential<real_t> exp1, exp2;

    explicit bimodal(
            const exponential<real_t> &exp1,
            const exponential<real_t> &exp2
    ) : exp1(exp1), exp2(exp2)
    {}

    real_t funval(const real_t lnr) const {
        return exp1(lnr) + exp2(lnr);
    }
};
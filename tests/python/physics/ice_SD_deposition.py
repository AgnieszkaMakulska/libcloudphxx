import sys

sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import lgrngn, common
import numpy as np
import csv

def run(T):

    backend = lgrngn.backend_t.serial
    opts_init = lgrngn.opts_init_t()
    opts = lgrngn.opts_t()

    kappa = .61
    soluble_fraction = 1.
    opts_init.dry_sizes = {(kappa, soluble_fraction): {0.01e-6: [1., 1]}}   # {r_dry: [concentration, number of SDs]}
    opts_init.n_sd_max = 1
    opts_init.dt = 0.1

    opts_init.ice_switch = True
    opts_init.init_ice = True
    opts_init.ice_r0 = 2e-6
    opts_init.coal_switch = False
    opts_init.sedi_switch = False
    opts.depo = True
    opts.cond = True
    opts.ice_nucl = False
    opts.coal = False
    opts.adve = False
    opts.sedi = False

    p = 101300.
    RH = 1.
    t_total = 600.

    # target values
    rv_target = RH * common.r_vs(T, p)
    th_target = T / common.exner(p)

    rv = np.array([rv_target])
    th = np.array([th_target])
    rhod = np.array([common.rhod(p, th[0], rv[0])])

    prtcls = lgrngn.factory(backend, opts_init)
    prtcls.init(th, rv, rhod)

    ice_a = np.array([])
    ice_c = np.array([])
    ice_rho = np.array([])

    for nt in range(int(t_total / opts_init.dt)):

        if(nt % 100 == 0):

            prtcls.diag_water()
            prtcls.diag_wet_mom(0)
            assert np.all(np.frombuffer(prtcls.outbuf()) == 0) # no liquid droplets

            prtcls.diag_ice()
            prtcls.diag_ice_a_mom(0)
            assert np.all(np.frombuffer(prtcls.outbuf()) == 1) # only one ice crystal

            ice_a = np.append(ice_a, np.asarray(prtcls.get_attr("ice_a")))
            ice_c = np.append(ice_c, np.asarray(prtcls.get_attr("ice_c")))
            ice_rho = np.append(ice_rho, np.asarray(prtcls.get_attr("ice_rho")))

        # reset the state to force constant T and RH
        th[0] = th_target
        rv[0] = rv_target

        prtcls.step_sync(opts, th, rv, rhod)
        prtcls.step_async(opts)

    return ice_a[-1], ice_c[-1], ice_rho[-1]

ice_a_array = np.array([])
ice_c_array= np.array([])
ice_rho_array = np.array([])
T_array = np.arange(243., 272., 0.1)
for T in T_array:
    ice_a, ice_c, ice_rho = run(T)
    ice_a_array = np.append(ice_a_array, ice_a)
    ice_c_array = np.append(ice_c_array, ice_c)
    ice_rho_array = np.append(ice_rho_array, ice_rho)

m_array = 4/3 * np.pi * ice_a_array * ice_a_array * ice_c_array * ice_rho_array
phi_array = ice_c_array / ice_a_array
m_ug = m_array * 1e9
a_um = ice_a_array * 1e6
c_um = ice_c_array * 1e6

with open("habit_vs_T_libcloud.csv", mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["T[K]", "m[ug]", "a[um]", "c[um]", "phi", "rho[kg/m3]"])
    for i in range(len(T_array)):
        writer.writerow([
            T_array[i],
            m_ug[i],
            a_um[i],
            c_um[i],
            phi_array[i],
            ice_rho_array[i]
        ])
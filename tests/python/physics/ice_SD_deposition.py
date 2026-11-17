import sys

sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import lgrngn, common
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)

backend = lgrngn.backend_t.serial
opts_init = lgrngn.opts_init_t()
opts = lgrngn.opts_t()

kappa = .61
soluble_fraction = 1.
opts_init.dry_sizes = {(kappa, soluble_fraction): {0.01e-6: [1., 1]}}   # {r_dry: [concentration, number of SDs]}
opts_init.n_sd_max = 10
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

p = 101000.
T = 265.
RH = 1.
t_total = 600

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
n_ice = np.array([])
n_liq = np.array([])

for nt in range(int(t_total / opts_init.dt)):
    # reset the state to force constant T and RH
    th[0] = th_target
    rv[0] = rv_target

    prtcls.step_sync(opts, th, rv, rhod)
    prtcls.step_async(opts)

    if(nt % 100 == 0):

        prtcls.diag_water()
        prtcls.diag_wet_mom(0)
        n_liq = np.append(n_liq, np.frombuffer(prtcls.outbuf()).copy())

        prtcls.diag_ice()
        prtcls.diag_ice_a_mom(0)
        n_ice = np.append(n_ice, np.frombuffer(prtcls.outbuf()).copy())

        ice_a = np.append(ice_a, np.asarray(prtcls.get_attr("ice_a")))
        ice_c = np.append(ice_c, np.asarray(prtcls.get_attr("ice_c")))
        ice_rho = np.append(ice_rho, np.asarray(prtcls.get_attr("ice_rho")))

assert np.all(n_liq == 0)
assert np.all(n_ice == 1)
print("ice_a [um]", ice_a * 1e6)
print("ice_c [um]", ice_c * 1e6)
print("ice_rho [kg/m3]", ice_rho)
time = np.arange(0, t_total, opts_init.dt * 100)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
ax1.plot(time / 60, ice_a * 1e6)
ax2.plot(time / 60, ice_c * 1e6)
ax3.plot(time / 60, ice_rho)
ax1.set_ylabel("a [um]")
ax2.set_ylabel("c [um")
ax3.set_ylabel("rho_i [kg/m3]")
ax3.set_xlabel("t [min]")
plt.tight_layout()
plt.savefig("ice_attributes.png")
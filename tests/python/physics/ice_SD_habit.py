import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15
plt.rcParams['lines.markersize'] = 4

df_scale = pd.read_csv('habit_vs_T_SCALE.csv')
T_scale = df_scale['T[K]'].to_numpy()
m_scale = df_scale['m[ug]'].to_numpy()
a_scale = df_scale['a[um]'].to_numpy()
c_scale = df_scale['c[um]'].to_numpy()
phi_scale = df_scale['phi'].to_numpy()
rho_i_scale = df_scale['rho[kg/m3]'].to_numpy()

df_libc = pd.read_csv('habit_vs_T_libcloud.csv')
T_libc = df_libc['T[K]'].to_numpy()
m_libc = df_libc['m[ug]'].to_numpy()
a_libc = df_libc['a[um]'].to_numpy()
c_libc = df_libc['c[um]'].to_numpy()
phi_libc = df_libc['phi'].to_numpy()
rho_i_libc = df_libc['rho[kg/m3]'].to_numpy()


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

ax1.scatter(T_scale, m_scale)
ax2.scatter(T_scale, a_scale, label="a")
ax2.scatter(T_scale, c_scale, label="c")
ax3.scatter(T_scale, rho_i_scale)

ax1.scatter(T_libc, m_libc)
ax2.scatter(T_libc, a_libc, label="a")
ax2.scatter(T_libc, c_libc, label="c")
ax3.scatter(T_libc, rho_i_libc)

ax2.legend(
    loc="upper left"
)
ax1.set_ylabel(r'm [$\mu$g]')
ax2.set_ylabel(r'a, c [$\mu$m]')
ax3.set_ylabel("apparent density [kg/m$^3$]")
ax3.set_xlabel("T [K]")
ax1.set_yscale("log")
ax2.set_yscale("log")
ax1.grid()
ax2.grid()
ax3.grid()
ax1.set_ylim([0.1, 20])
ax2.set_ylim([1, 2e3])
ax3.set_ylim([250, 950])
plt.tight_layout()
plt.savefig("ice_habits.png")



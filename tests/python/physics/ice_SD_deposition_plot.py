import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15
plt.rcParams['lines.markersize'] = 3

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

# Comparison with SCALE-SDM

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

ax1.scatter(T_libc, m_libc, label='UWLCM', color = 'black')
ax2.scatter(T_libc, a_libc, color = 'tab:blue', label="a-axis UWLCM")
ax2.scatter(T_libc, c_libc, color = 'tab:orange', label="c-axis UWLCM")
ax3.scatter(T_libc, rho_i_libc, label='UWLCM', color = 'black')

ax1.scatter(T_scale, m_scale, label = 'SCALE-SDM', color = 'red', marker = 'o')
ax2.scatter(T_scale, a_scale, color = 'tab:cyan', label="a-axis SCALE-SDM", marker = 'o')
ax2.scatter(T_scale, c_scale, color = 'tab:red', label="c-axis SCALE-SDM", marker = 'o')
ax3.scatter(T_scale, rho_i_scale, label = 'SCALE-SDM', color = 'red', marker = 'o')

ax1.legend()
ax2.legend(
    ncols = 2,
    loc="lower center"
)
ax3.legend()
ax1.set_ylabel(r'm [$\mu$g]')
ax2.set_ylabel(r'a, c [$\mu$m]')
ax3.set_ylabel(r'$\rho_{app}$ [kg/m$^3$]')
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


# Relative difference between models

df_sub = df_libc[df_libc['T[K]'].round(1) % 1 == 0]
T_libc = df_sub['T[K]'].to_numpy()
m_libc = df_sub['m[ug]'].to_numpy()
a_libc = df_sub['a[um]'].to_numpy()
c_libc = df_sub['c[um]'].to_numpy()
phi_libc = df_sub['phi'].to_numpy()
rho_i_libc = df_sub['rho[kg/m3]'].to_numpy()

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 12), sharex=True)
ax1.scatter(T_libc, (m_libc - m_scale) / m_scale, color = 'black')
ax2.scatter(T_libc, (a_libc - a_scale) / a_scale, color = 'black')
ax3.scatter(T_libc, (c_libc - c_scale) / c_scale, color = 'black')
ax4.scatter(T_libc, (rho_i_libc - rho_i_scale) / rho_i_scale, color = 'black')
ax1.set_ylabel(r'$\frac{m_{UWLCM} - m_{SCALE}}{m_{SCALE}}$')
ax2.set_ylabel(r'$\frac{a_{UWLCM} - a_{SCALE}}{a_{SCALE}}$')
ax3.set_ylabel(r'$\frac{c_{UWLCM} - c_{SCALE}}{c_{SCALE}}$')
ax4.set_ylabel(r'$\frac{\rho_{i,UWLCM} - \rho_{i,SCALE}}{\rho_{i,SCALE}}$')
ax4.set_xlabel("T [K]")
for ax in fig.axes:
    ax.yaxis.label.set_size(20)
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
plt.tight_layout()
plt.savefig("ice_habits_difference.png")
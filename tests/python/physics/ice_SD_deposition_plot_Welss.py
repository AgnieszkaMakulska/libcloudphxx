import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15
plt.rcParams['lines.markersize'] = 2

# UWLCM (libcloudph++) runs: (file, label, color, c-axis color)
libc_runs = [
    ('habit_vs_T_libcloud_Shima_2020.csv',                           'Shima 2020',                                                  'black',      'tab:orange'),
    #('habit_vs_T_libcloud_Welss_IGF.csv',                           'Welss IGF',                                                   'tab:blue',   'tab:pink'),
    ('habit_vs_T_libcloud_Welss_IGF_branching.csv',                  'Welss IGF2 + branching above 200 um',                         'tab:green',  'tab:red'),
    #('habit_vs_T_libcloud_Welss_IGF_branching_ventilation.csv',     'Welss IGF + branching + vent.',                               'tab:orange', 'tab:brown'),
    ('habit_vs_T_libcloud_Welss_IGF_branching_ventilation_r_eq.csv', 'Welss IGF2 + branching above 200 um + modified ventilation ', 'tab:purple', 'tab:cyan'),
]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

for fname, label, color, color_c in libc_runs:
    df_libc = pd.read_csv(fname)
    T_libc = df_libc['T[K]'].to_numpy()
    m_libc = df_libc['m[ug]'].to_numpy()
    a_libc = df_libc['a[um]'].to_numpy()
    c_libc = df_libc['c[um]'].to_numpy()
    phi_libc = df_libc['phi'].to_numpy()
    rho_i_libc = df_libc['rho[kg/m3]'].to_numpy()

    ax1.scatter(T_libc, m_libc, label=label, color=color)
    ax2.scatter(T_libc, a_libc, color=color, marker='o', label="a-axis " + label)
    ax2.scatter(T_libc, c_libc, color=color_c, marker='x', label="c-axis " + label)
    ax3.scatter(T_libc, rho_i_libc, label=label, color=color)

ax1.legend(fontsize=9)
ax2.legend(
    ncols = 2,
    loc="lower center",
    fontsize=7
)
ax3.legend(fontsize=9)
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
plt.savefig("ice_habits_Welss.png", dpi=300)
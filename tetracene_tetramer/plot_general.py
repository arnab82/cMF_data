import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors

#copy the energy values from the txt file

BFGS = []
CG = []
diis_gd = []
diis_wop = []
diis_wp = []
newton_wop = []
newton_wp = []
# print(BFGS)

# Generate sample data for different methods
iteration_bfgs = np.arange(2, len(BFGS) + 2)
print(iteration_bfgs)
iteration_cg = np.arange(2, len(CG) + 2)
# iteration_gd = np.arange(1, len(GD) + 1)
iteration_diis_gd = np.arange(2, len(diis_gd) + 2)
iteration_diis_wop = np.arange(2, len(diis_wop) + 2)
iteration_diis_wp = np.arange(2, len(diis_wp) + 2)
iteration_newton_wop= np.arange(2, len(newton_wop) + 2)
iteration_newton_wp= np.arange(2, len(newton_wp) + 2) 
blue = '#3E6D9C'
orange = '#FD841F'
red_orange = '#E14D2A'
dark_blue = '#001253'
cb = ['#000000']
cb.extend([i for i in plt.rcParams['axes.prop_cycle'].by_key()['color']])
# print(CG)
# Create a plot
# print(np.array(BFGS) - BFGS[-1])
plt.grid(False)
plt.plot(iteration_bfgs, np.log10(np.array(BFGS) - BFGS[-1]), label='BFGS', marker='o', markersize=4, color=cb[0])
plt.plot(iteration_cg, np.log10(np.array(CG) - CG[-1]), label='CG', marker='s', markersize=4, color=cb[1])
plt.plot(iteration_diis_gd, np.log10(np.array(diis_gd) - diis_gd[-1]), label='DIIS_GD', marker='+', markersize=4, color=cb[2])
plt.plot(iteration_diis_wop, np.log10(np.array(diis_wop) - diis_wop[-1]), label='DIIS_WOP', marker='^', markersize=4, color=cb[3])
plt.plot(iteration_diis_wp, np.log10(np.array(diis_wp) - diis_wp[-1]), label='DIIS_WP', marker='p', markersize=4, color=cb[4])
plt.plot(iteration_newton_wop, np.log10(np.array(newton_wop) - newton_wop[-1]), label='NEWTON_WOP', marker='x', markersize=4, color=cb[5])
plt.plot(iteration_newton_wp, np.log10(np.array(newton_wp) - newton_wp[-1]), label='NEWTON_WP', marker='D', markersize=4, color=cb[6])
# plt.xscale('symlog', linthresh=5.0)

# plt.plot(iteration_bfgs, np.log10(np.array(g_bfgs) - g_bfgs[-1]), label='BFGS', marker='o')
# plt.plot(iteration_cg, np.log10(np.array(g_cg) - g_cg[-1]), label='CG', marker='s')
# plt.plot(iteration_diis_gd, np.log10(np.array(g_diis_gd) - g_diis_gd[-1]), label='DIIS_GD', marker='+')
# plt.plot(iteration_diis_wop, np.log10(np.array(g_diis_wop) - g_diis_wop[-1]), label='DIIS_WOP', marker='x')
# plt.plot(iteration_diis_wp, np.log10(np.array(g_diis_wp) - g_diis_wp[-1]), label='DIIS_WP', marker='p')
# plt.plot(iteration_newton_wop, np.log10(np.array(g_newton_wop) - g_newton_wop[-1]), label='NEWTON_WOP', marker='^')
# plt.plot(iteration_newton_wp, np.log10(np.array(g_newton_wp) - g_newton_wp[-1]), label='NEWTON_WP', marker='D')
# Customize the plot
# plt.title('Gradient vs. Iteration Number for Different Methods')
# plt.title('Energy vs. Iteration Number for Different Methods')
# plt.ylabel('Gradient')
plt.ylabel('Log$_{10}$ (Energy \N{MINUS SIGN} Converged Energy) ', fontsize='12')
plt.xlabel('Iteration Number', fontsize='12')
plt.title("Comparison between cMF methods for $[$Fe$_2$OCl$_6]^{2-}$") 
plt.legend()
fig = plt.gcf()
fig.set_size_inches(5, 5)
# fig.savefig("fe_energy_logscale.pdf", dpi=600, bbox_inches='tight')

plt.show()

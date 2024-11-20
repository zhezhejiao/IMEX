

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


plt.style.use('classic')
plt.rc('font', family='Times New Roman', size=13)

#Ensure the file path is correctly modified!
results_folder = './results/rIMEX'

data = pd.read_csv(results_folder + ".txt", delim_whitespace=True)

plt.figure(figsize=(8, 6))

tau_ref_2 = data.iloc[1, 2]
error_ref_2 = data.iloc[1, 3]

C_2 = error_ref_2 / (tau_ref_2**2)

x_vals = np.logspace(np.log10(1e-5), np.log10(1e0), 500)

convergence_curve_2 = C_2 * x_vals**2

tau_ref_1 = data.iloc[-3, 2]
error_ref_1 = data.iloc[-3, 3]

C_1 = error_ref_1 / tau_ref_1

convergence_curve_1 = C_1 * x_vals

plt.plot(x_vals, convergence_curve_2, label='order 2', linestyle='--', linewidth=2, color='black')

plt.plot(x_vals, convergence_curve_1, label='order 1', linestyle=':', linewidth=2, color='gray')

plt.plot(data.iloc[:, 2], data.iloc[:, 3], label='rIMEX error',
         marker='x', markersize=8, markerfacecolor='none', markeredgewidth=2, markeredgecolor='maroon',
         linestyle='-', linewidth=2, color='red', alpha=0.9)

plt.legend(loc='best')

plt.xlabel('$\\tau$', fontsize=14, labelpad=-4)
plt.ylabel('Error', fontsize=14)

plt.xscale('log')
plt.yscale('log')

plt.xlim(1e-4, 1e-1)
plt.ylim(1e-6, 1e0)

plt.tick_params(which='both', direction='in', length=6, width=1, colors='black',
               grid_color='gray', grid_alpha=0.5, top=True, right=True, labelsize=10)
plt.tick_params(which='minor', axis='y', left=False, right=False)

plt.show()

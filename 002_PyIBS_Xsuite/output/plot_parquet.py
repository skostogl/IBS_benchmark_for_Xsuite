import numpy as np
import matplotlib.pylab as plt
import pandas as pd

kinetic    = pd.read_parquet("xsuite_kinetic.parquet")
#simple     = pd.read_parquet("xsuite_simple.parquet")
analytical = pd.read_parquet("xsuite_analytical.parquet")

f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (16,5))

plt.sca(ax1)
plt.plot(kinetic['eps_x'].values*1e6, c='b', label='kinetic')
#plt.plot(simple['eps_x'].values*1e6, c='g', label='simple kick')
plt.plot(analytical['eps_x'].values*1e6, c='k', label='analytical')
plt.legend(fontsize=12)

plt.sca(ax2)
plt.plot(kinetic['eps_y'].values*1e6, c='b')
#plt.plot(simple['eps_y'].values*1e6, c='g')
plt.plot(analytical['eps_y'].values*1e6, c='k')

plt.sca(ax3)
plt.plot(kinetic['sig_delta'].values*1e3, c='b')
#plt.plot(simple['sig_delta'].values*1e3, c='g')
plt.plot(analytical['sig_delta'].values*1e3, c='k')

plt.sca(ax4)
plt.plot(kinetic['bl'].values*1e2, c='b')
#plt.plot(simple['bl'].values*1e2, c='g')
plt.plot(analytical['bl'].values*1e2, c='k')

ax1.set_ylabel(r'$\varepsilon_x$ [um]')
ax1.set_xlabel('Turns')

ax2.set_ylabel(r'$\varepsilon_y$ [um]')
ax2.set_xlabel('Turns')

ax3.set_ylabel(r'$\sigma_{\delta}$ [$10^{-3}$]')
ax3.set_xlabel('Turns')

ax4.set_ylabel(r'bl [cm]')
ax4.set_xlabel('Turns')

plt.tight_layout()

plt.savefig('comparison_parquet.png', dpi = 400)
plt.show()

import scipy.constants
from scipy.constants import e as qe
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_python = pd.read_parquet("IBS_output_python.parquet") 
df_madx = pd.read_parquet("IBS_output_madx.parquet")

tt_p = df_python["tt"]
tt_m = df_madx["tt"]

fig, ax = plt.subplots(ncols=3, figsize=(10,6))
plt.sca(ax[0])
plt.ylabel("emitx (um)")

plt.plot(tt_p,df_python["exin"].values*1e6, marker='o', c='b', label='Python IBS', lw=4, ms=10)
plt.plot(tt_m,df_madx["exin"].values*1e6, marker='o', c='r', label='MADX IBS', lw=3, ms=5)

plt.legend()

plt.sca(ax[1])
plt.ylabel("emity (um)")
plt.plot(tt_p,df_python["eyin"].values*1e6, marker='o', c='b', lw=4, ms=10)
plt.plot(tt_m,df_madx["eyin"].values*1e6, marker='o', c='r', lw=3, ms=5)


plt.sca(ax[2])
plt.ylabel("bunch length (cm)")
plt.plot(tt_p, 1e2*df_python["bl_ns"].values*1e-9*scipy.constants.speed_of_light/4., marker='o', c='b', lw=4, ms=10)
plt.plot(tt_m, 1e2*df_madx["bl_ns"].values*1e-9*scipy.constants.speed_of_light/4., marker='o', c='r', lw=3, ms=5)
fig.text(0.5, 0.0, 'Time (sec.)', ha='center')

fig.tight_layout()
fig.savefig("example.png")
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import plotly.express as px
import seaborn as sns

sns.set('paper')

## old data sets
# path_ricc2='/Users/carlton/Work/ADCHAM_ClusterIn/output/RICC2/'
# path_dlpno='/Users/carlton/Work/ADCHAM_ClusterIn/output/DLPNO/'
# path_old='/Users/carlton/Work/ADCHAM_ClusterIn/output/old_ACDC_1/'

path_output='/Users/carlton/Work/Clusterin_multiple_chemistries/ClusterIn/ClusterIn_plugin/output/'


Jnucl_N=np.genfromtxt(path_output+'Jnucl_N.dat')
Jnucl_D=np.genfromtxt(path_output+'Jnucl_D.dat')
Jnucl_total=np.genfromtxt(path_output+'Jnucl_total.dat')

# DLPNO=np.genfromtxt('/Users/carlton/Work/ADCHAM_ClusterIn/ClusterIn_plugin/output/Jnucl_traditional.dat')
# RICC2=np.genfromtxt('/Users/carlton/Work/ADCHAM_ClusterIn/ClusterIn_plugin_old_acdc_data/output/Jnucl_traditional.dat')

fig=plt.subplots(nrows=1,ncols=1,figsize=(6,4))
plt.subplot(111)

plt.plot(Jnucl_N, color='r', linestyle='-', linewidth=2, label='NH3')
plt.plot(Jnucl_D, color='k', linestyle='-', linewidth=2, label='DMA')
plt.plot(Jnucl_total, color='b', linestyle='--', linewidth=2, label='Total')
plt.yscale('log')
plt.grid('True')
plt.xlabel('Simulation Time (min)')
plt.title(r'$J for different chemistries$')
plt.ylabel('J')
plt.legend()
# plt.ylim([1e-10,max(Jnucl_ricc2_acdc)])
plt.show()
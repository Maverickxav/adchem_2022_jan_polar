import numpy as np
import matplotlib.pyplot as plt

path='/Users/carlton/Work/ACDC_versions/Clusterin_multiple_chemistries/ClusterIn/ClusterIn_msa/output/'

jnucl_D=np.genfromtxt(path+'Jnucl_D.dat')
jnucl_I=np.genfromtxt(path+'Jnucl_I.dat')
jnucl_M=np.genfromtxt(path+'Jnucl_M.dat')
jnucl_N=np.genfromtxt(path+'Jnucl_N.dat')
jnucl_tot=np.genfromtxt(path+'Jnucl_total.dat')
#

fig = plt.subplots(1,1,figsize=(8,6))

plt.subplot(111)
plt.plot(np.arange(0,len(jnucl_D)), jnucl_D,c='indianred', label='AD')
plt.plot(np.arange(0,len(jnucl_D)), jnucl_I,c='dodgerblue', label='IoIi')
plt.plot(np.arange(0,len(jnucl_D)), jnucl_M,c='magenta', label='AMsD')
plt.plot(np.arange(0,len(jnucl_D)), jnucl_N,c='slategrey', label='AN')
# plt.plot(np.arange(0,len(jnucl_D)), jnucl_tot,c='k', label='Total')
plt.legend()
plt.ylabel('Jnucl')
plt.yscale('log')
plt.grid(True)
plt.xlabel('no. of iterations')
plt.show()

# plt.subplot(212)
# plt.plot(np.arange(0,len(jnucl_ad)), cout_ad,c='indianred', label='AD')
# plt.plot(np.arange(0,len(jnucl_ad)), cout_amd,c='dodgerblue', label='AMD')
# plt.legend()
# plt.ylabel('cout')
# plt.yscale('log')
# plt.grid(True)
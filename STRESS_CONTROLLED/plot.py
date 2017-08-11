# Librairies necessaires
import matplotlib        
import matplotlib.pyplot as plt
import numpy 			 as np
import seaborn           as sns
import pylab 			 as p
import sys



pp_effects = True


# Layers 1 and 2
m1 = 0.5736
m2 = 0.4067


# Loading ppts
tstep = 0.002

# Files
file_sigeps   = 'outputfiles/stressstrainxz'
file_backbone = 'outputfiles/SOILHYPER'
file_gref     = 'outputfiles/gref'
file_Gmod     = 'outputfiles/Gmodulus'
file_S0       = 'outputfiles/S0'
file_w        = 'outputfiles/w'
file_deps     = 'outputfiles/depsxz'
file_S        = 'outputfiles/deviatoriceffective2'


# Reading files
gamma    = np.genfromtxt(file_sigeps,usecols=0)
sigma    = np.genfromtxt(file_sigeps,usecols=1)

bb_gamma = np.genfromtxt(file_backbone,usecols=0)
bb_sigma = np.genfromtxt(file_backbone,usecols=1)

# Respecting loading frequency in main.f90
time     = np.linspace(0.0,len(sigma)*tstep/4.0, num=len(sigma))
print 'Total time   ', len(sigma)*tstep
print '1 cycle      ', len(sigma)*tstep/4.0


gref     = np.genfromtxt(file_gref,usecols=1)
deps     = np.genfromtxt(file_deps,usecols=0)



###
print 
print 'Maximum strain [%]    : ', max(abs(gamma))*1.1
print 'Applied Stress [kPa]  : ', max(sigma)

print 'Initial Reference strain : ', gref[0]





###
if pp_effects:
	Gmod     = np.genfromtxt(file_Gmod,usecols=1)
	print 'Minimum G modulus:'
	print min(Gmod)

	S0       = np.genfromtxt(file_S0,usecols=1)
	w        = np.genfromtxt(file_w, usecols=0)
	S        = np.genfromtxt(file_S, usecols=0)
	r        = np.genfromtxt(file_S, usecols=1)


	# Failure and phase transformation lines
	lfailure = np.zeros((15))
	lphase   = np.zeros((15))
	array    = np.zeros((15))

	for i in np.arange(0,1.5,0.1):
	    lfailure[i] = i* m1
	    lphase  [i] = i* m2
	    array   [i] = i
###





###
print 'Plotting'
# fig = p.figure(figsize=(18,6))
fig = p.figure(figsize=(11.69,8.27))

p.subplots_adjust(hspace=0.35)
sns.set_style('whitegrid')


# Stress-strain curve
ax = fig.add_subplot(231)

xlimit = max(abs(gamma))*1.1
ylimit = max(abs(sigma))*1.1

ax.set_xlim([-xlimit,xlimit])
ax.set_ylim([-ylimit,ylimit])

ax.plot(gamma,sigma, color='red')
ax.plot(bb_gamma,bb_sigma, color='black')
ax.plot(-bb_gamma,-bb_sigma, color='black')

ax.set_xlabel('Shear strain [%]', fontsize=14)
ax.set_ylabel('Shear stress [kPa]', fontsize=14)

plt.xticks(fontsize=13)
plt.yticks(fontsize=13)


# Stress vs time
ax = fig.add_subplot(232)
# plt.title('Layer 1 - CSR=0.15'+ '\n', fontsize=16)
plt.title('Layer 2 - CSR=0.20'+ '\n', fontsize=16)


ax.plot(time,sigma, color='black')
ax.set_xlabel('Number of loading cycles', fontsize=14)
ax.set_ylabel('Stress [kPa]', fontsize=14)

plt.xticks(fontsize=13)
plt.yticks(fontsize=13)


#
# ax = fig.add_subplot(233)
# #ax.plot(time,gref, color='black', marker='o')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Reference strain')


#
ax = fig.add_subplot(234)
# ax.set_xscale('log')
# ax.plot(time,deps, color='black', marker='o')
ax.plot(time,gamma/2.0, color='black')

ax.set_xlabel('Number of loading cycles', fontsize=14)
ax.set_ylabel('Axial strain [%]', fontsize=14)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)





if pp_effects:

	#
	ax = fig.add_subplot(233)
	ax.set_xlim([0,1.2])
	ax.plot(S,r, color='black')
  	ax.plot(array, lfailure, color = 'gray')
  	ax.plot(array, lphase, color = 'gray')
  	ax.plot(array, lphase*0.67,   color = 'gray', linestyle='--')


	ax.set_xlabel('Normalized effective stress', fontsize=14)
	ax.set_ylabel('Normalized deviatoric stress', fontsize=14)
	plt.xticks(fontsize=13)
	plt.yticks(fontsize=13)

	#
	ax = fig.add_subplot(235)
	ax.plot(time,S0, color='black')
	ax.set_xlabel('Number of loading cycles', fontsize=14)
	ax.set_ylabel('S0', fontsize=14)
	plt.xticks(fontsize=13)
	plt.yticks(fontsize=13)

	#
	ax = fig.add_subplot(236)
	ax.plot(time,Gmod, color='black')
	ax.set_xlabel('Number of loading cycles', fontsize=14)
	ax.set_ylabel('Shear modulus [MPa]', fontsize=14)
	plt.xticks(fontsize=13)
	plt.yticks(fontsize=13)
#


plt.tight_layout()
# fig.savefig('Layer1_csr0_15.png',dpi=300)
fig.savefig('Layer2_csr0_20.png',dpi=300)

plt.show()

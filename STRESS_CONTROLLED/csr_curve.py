# SINGLE PLOTS
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import pylab as p
import numpy as np
import seaborn as sns
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

# Files to compare with

# Layer 1
n1 = [21,7,3]
csr1   = [0.15,0.20,0.28]

# Layer 2
n2 = [12,3.5,2]
csr2   = [0.20,0.30,0.40]






#Fig style
print 'Plotting'
fig = p.figure(figsize=(8,6))
p.subplots_adjust(hspace=0.35)
sns.set_style('whitegrid')


#
ax1 = fig.add_subplot(111)

ax1.set_xscale('log')
ax1.set_xlim([1,100])
ax1.set_ylim([0,0.45])

ax1.plot(n1,csr1, marker = 'o', label = 'Layer 1', color='red')
ax1.plot(n2,csr2, marker = 'o', label = 'Layer 2', color='black')
ax1.legend(prop={'size':13})

plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid(True, which='both')

ax1.set_xlabel('Number of loading cycles', fontsize=18)
ax1.set_ylabel('Cyclic stress ratio', fontsize=18)




p.savefig("csr.png",dpi=300)
plt.show()
#
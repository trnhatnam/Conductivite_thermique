import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

#os.chdir("C:/Users/trinh/Desktop/Project_elec")

meijerg = pd.read_csv('meijerg_only.csv')
simpson = pd.read_csv('simpson.csv')

diff_re = np.abs(meijerg['Re'] - simpson['Re'])
diff_im = np.abs(meijerg['Im'] - simpson['Im'])
freq = meijerg['Thermal_freq']

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.set_title('Différence entre les deux méthodes sur la température moyenne (silicium 5e-6)')
ax1.semilogx(freq,diff_re, color="blue",label="réel")
ax1.semilogx(freq,diff_im, color="orange",label="imaginaire")
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_xlabel('Thermal_freq')
ax1.set_ylabel('Différence de température')
ax1.legend()


ax2 = fig.add_subplot(122)
ax2.set_title('Superposition Meijerg-Simpson (silicium 5e-6)')
ax2.semilogx(freq,meijerg['Re'], color="blue", label="réel meij")
ax2.semilogx(freq,simpson['Re'], color="orange",label="réel simp")
ax2.semilogx(freq,meijerg['Im'], color="blue", linestyle="--",label="imag meij")
ax2.semilogx(freq,simpson['Im'], color="orange", linestyle="--", label="imag simp")
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_xlabel('Thermal_freq')
ax2.set_ylabel('Température moyenne')
ax2.legend()

plt.show()

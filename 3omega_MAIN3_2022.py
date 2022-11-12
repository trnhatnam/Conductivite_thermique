from mpmath import meijerg, j
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt


# pip install matplotlib
# pip install numpy

# x = np.arange(0,4*np.pi,0.1)   # start,stop,step
# y = np.sin(x)
# z = np.cos(x)
# plt.plot(x,y,x,z)
# plt.xlabel('x values from 0 to 4pi')  # string must be enclosed with quotes '  '
# plt.ylabel('sin(x) and cos(x)')
# plt.title('Plot of sin and cos from 0 to 4pi')
# plt.legend(['sin(x)', 'cos(x)'])      # legend entries as seperate strings in a list
# plt.show()


# parametres de calcul 3 Omega
#materiau
D = 1.3e-7      # diffusivité  (intrinsèque au materiau etudié)
kappa = 0.3     # conductivité thermique (intrinsèque au materiau etudié)
#Resistor
TCR = 0.00253   # coefficient de temperature resistor (Au) - (intrinsèque au materiau du resistor)
L = 0.002       # longueur du resistor
V0 = 1          # pic de tension fondamental
R0 = 80         # R resistor à temperature ambiante
P = xxx         # pissance par unité de longueur (w/m)

gamma = 0.5772  #constante d`Euler

#variables
ts = np.array([350e-6])                         # epaisseur du substrat
bh = np.array([5e-6, 10e-6, 15e-6, 20e-6])      #demi largeur du resistor
frequence = np.arange(1, 100000, 1000)          #domaine d`étude frequentiel


#Etude asymptotique

#Resolution numerique eq 4 ref2 - methode des fonctions hypergeometriques MeijerG

#Resolution numerique eq 4 ref2 - integration numerique methode de Simpson

#Representation graphique


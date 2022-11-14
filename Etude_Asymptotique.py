from mpmath import meijerg, j
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt



# parametres de calcul 3 Omega
#materiau
D = 1.3e-7      # diffusivité  (intrinsèque au materiau etudié)
kappa = 0.3     # conductivité thermique (intrinsèque au materiau etudié)
#Resistor
TCR = 0.00253   # coefficient de temperature resistor (Au) - (intrinsèque au materiau du resistor)
L = 0.002       # longueur du resistor
V0 = 1          # pic de tension fondamental
R0 = 80         # R resistor à temperature ambiante
P = 1         # pissance par unité de longueur (w/m)

gamma = 0.5772  #constante d`Euler

#variables
ts = np.array([350e-6])                         # epaisseur du substrat
bh = np.array([5e-6, 10e-6, 15e-6, 20e-6])      #demi largeur du resistor
frequence = np.arange(1, 100000, 1000)          #domaine d`étude frequentiel
epsilon = 0.923  #Donnée dans la référence diapo 33 ref1


#Etude asymptotique
for b in bh :
    Tac = -P/(2*np.pi*kappa)*(np.log(b**2/D)+np.log(4*np.pi*frequence)-2*epsilon)-1j*P/(4*kappa)
    Module = np.abs(Tac) #Module 
    Argument = np.angle(Tac) #Argument 

    #In-phase 
    Tac_inphase=Module*np.cos(Argument)
    plt.plot(np.log(4*np.pi*frequence),Tac_inphase,label="b="+str(b))


#Out-phase
Tac_outphase =Module*np.sin(Argument)

plt.plot(np.log(4*np.pi*frequence),Tac_outphase,"--")
plt.xlabel("ln(2w)")
plt.ylabel("ΔTac")
plt.legend(loc="upper right")
plt.title("Asymptote à haute fréquence de  ΔTac en fonction de ln(2w)")
plt.show()


#Resolution numerique eq 4 ref2 - methode des fonctions hypergeometriques MeijerG

#Resolution numerique eq 4 ref2 - integration numerique methode de Simpson

#Representation graphique


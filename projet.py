from mpmath import meijerg, j
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys # ajouté


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


#Etude asymptotique

# Fonctions de conversion
def no_div_by_0(arr):
    # Remplace les 0 par un epsilon dans l'array
    cp = arr.copy()
    for (i,donnee) in enumerate(arr):
        if abs(donnee) < sys.float_info.epsilon:
            cp[i] = sys.float_info.epsilon
    return cp

def lamb_2_ln2w(lamb):
    # Pronfdeur de pénétration -> ln(2w)
    lamb = no_div_by_0(lamb)
    return np.log(D/(lamb*lamb))

def ln2w_2_lamb(ln2w):
    # ln(2w) -> Profondeur de pénétration
    ln2w = no_div_by_0(ln2w)
    return np.sqrt(D/np.exp(ln2w))

# Plot pour la basse fréquence
fig, ax1 = plt.subplots()
pulsation = 2*np.pi*frequence
ln2w = np.log(2*pulsation)
eps = 0.923

for b in bh:
    # Calcul du deltaAc complexe
    deltaAc = (-P/(2*np.pi*kappa)) * (np.log(b*b/D) + ln2w - 2*eps) - 1j*P/(4*kappa) # basse fréquence (1.51 ref 2)
    #deltaAc2 = (P/(2*b*kappa*np.sqrt(np.exp(ln2w)/D)))*np.exp(-1j*np.pi/4) # haute fréquence (1.56 ref 2)
    
    # Module et arg
    moduleDeltaAc = np.abs(deltaAc)
    argDeltaAc = np.angle(deltaAc)
    #moduleDeltaAc2 = np.abs(deltaAc2)
    #argDeltaAc2 = np.angle(deltaAc2)
    
    # Plot
    p = plt.plot(ln2w, moduleDeltaAc*np.cos(argDeltaAc), '-', label="hf en phase b = " + str(b)) # en phase basse fréquence
    #plt.plot(ln2w, moduleDeltaAc2*np.cos(argDeltaAc2), ':', color=p[-1].get_color()) # en phase haute fréquence
    #plt.plot(ln2w, moduleDeltaAc2*np.sin(argDeltaAc2), '-.') # quadrature haute fréquence


# Quadrature basse fréquence (indépendant de b)
ax1.plot(ln2w, moduleDeltaAc*np.sin(argDeltaAc), '--', label="en quadrature")

# Profondeur de pénétration
ax2 = ax1.secondary_xaxis('top', functions=(ln2w_2_lamb, lamb_2_ln2w)) 
ax2.set_xlabel('Profondeur de pénétration (en m)')
ax2.invert_xaxis() # pas linéairement indépendant à ln2w
ln2w_ticks = ax1.get_xticks()
lamb_ticks = ln2w_2_lamb(ln2w_ticks) # coincider ln(2w) avec la profondeur de pénétration
ax2.set_xticks(lamb_ticks)
ax2.ticklabel_format(axis="x", style="sci", scilimits=(0,0)) # écriture scientifique

# Estéthique
ax1.set_title('Asymptotes de la température en fonction du log de la pulsation')
ax1.set_xlabel('ln(2w)')
ax1.set_ylabel('deltaAC (en °C)')
plt.legend(loc="upper right")
plt.show()

#Resolution numerique eq 4 ref2 - methode des fonctions hypergeometriques MeijerG

#Resolution numerique eq 4 ref2 - integration numerique methode de Simpson

#Representation graphique

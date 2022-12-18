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
k = 0.3     # conductivité thermique (intrinsèque au materiau etudié)
#Resistor
TCR = 0.00253               # coefficient de temperature resistor (Au) - (intrinsèque au materiau du resistor)
L = 0.002                   # longueur du resistor
V0 = 1                      # pic de tension fondamental
R0 = 80                     # R resistor à temperature ambiante
P = V0**2/(2*R0)            # pissance par unité de longueur (w/m)

gamma = 0.5772  #constante d`Euler

#variables
ts = np.array([350e-6, 400e-6])                         # epaisseur du substrat
bh = np.array([5e-6, 10e-6, 20e-6])      #demi largeur du resistor
frequence = np.arange(1, 1000, 200)          #domaine d`étude frequentiel

# Cartesian product of input variables
idx = pd.MultiIndex.from_product([bh, ts, frequence], names=["bh", "ts", "frequence"])
# print(idx)
# Reset the index of the DataFrame, and use the default one instead. If the DataFrame has a MultiIndex, this method can remove one or more levels.
df = pd.DataFrame(index=idx).reset_index()
# print(df)

omega = (df["bh"].values**2 * df["frequence"].values) * (2 * math.pi / D)     # Omega is vectorized naturally.
thermal_freq = (2 * df['frequence'])                                                # Thermal frequency is 2nd harmonic
T_depth = (np.sqrt(2*D / (2 * (math.pi) * thermal_freq)))/1e-6                      # Thermal penetration depth function of omega themal not electrical - Cahill
                                                                                    # sqrt(2D/wt)
# print(T_depth)

def f_u(omega_elem):
     
    asympt = (P / (k*L*math.pi)) * (-(1 / 2) * np.log(omega_elem) + 3 / 2 - gamma - j * ((math.pi) / 4))

    #ajouter fonctions exacxtes calculée en fonction de MeijerG ou Simpson
    # meijerg
    val1 = (-j*P / (4*L*k*math.pi * omega_elem)) * meijerg([[1, 3 / 2], []], [[1, 1], [0.5, 0]], j * omega_elem)  #solution approximée via fnction MeijerG on recupere reel et imaginaire

    #faire les calculs ici : amplitude , phase, V3omega
    amplitude = math.sqrt(np.real(val1) ** 2 + np.imag(val1) ** 2)
    phase = math.degrees(math.atan(np.imag(val1) / np.real(val1)))
    V3omega_asympt = 0.5 * V0 * TCR * asympt                                        # calculate thrid harmonic from DT
    V3omega = 0.5 * V0 * TCR * val1            

    return val1, asympt, amplitude, phase, V3omega, V3omega_asympt   #ajouter fle retour des calculs


#return a tuple of array. Remember to assign two otypes.
f_u_vec = np.vectorize(f_u, otypes=[np.complex128,      # MeijerG is complex number ou Simpson
                                    np.complex128,      # asympt is complex number
                                    np.ndarray,         # amplitude is array
                                    np.ndarray,         # phase is array
                                    np.complex128,      # V3omega_asympt is complex number
                                    np.complex128]      # V3omega is complex number issu de MeijerG ou Simpson
                                    )

tup = f_u_vec(omega) # tuple of arrays: (val1, asympt, amplitude, phase, V3omega_asympt, V3omega)
meij, asympt, amplitude, phase, V3omega, V3omega_asympt  = tup

### Fait par Nhat-nam le 18/12
# temp vs freq_elec
for choixB in range(0,1): # remplacer par range(bh.size) pour parcourir tous les b (attention c'est illisible)
    # extraction des valeurs qui vont nous intéresser
    choixIdx = df.index[(df["bh"] == bh[choixB])&(df["ts"]==ts[0])].to_numpy() # retourne les index des lignes qui vérifient b = bh[choixTs] et ts = ts[choixTs]
    omegaEtude = np.take(omega, choixIdx)
    deltaT_in_M = np.take(np.real(meij), choixIdx)
    deltaT_out_M = np.take(np.imag(meij), choixIdx)
    
    frequence_elec = (omegaEtude*D)/(4*np.pi*bh[choixB]**2) # on revient à la fréquence électrique à partir de omega
    plt.semilogx(frequence_elec, deltaT_in_M)

plt.title("delta_T en fonction de la frequence électrique")
plt.xlabel("Frequence électrique")
plt.ylabel("Température moyenne (en °C)")
plt.show()

# temp vs frequence thermique
for choixB in range(0,1): # remplacer par range(bh.size) pour parcourir tous les b (attention c'est illisible)
    # extraction des valeurs qui vont nous intéresser
    choixIdx = df.index[(df["bh"] == bh[choixB])&(df["ts"]==ts[0])].to_numpy() # retourne les index des lignes qui vérifient b = bh[choixTs] et ts = ts[choixTs]
    omegaEtude = np.take(omega, choixIdx)
    deltaT_in_M = np.take(np.real(meij), choixIdx)
    deltaT_out_M = np.take(np.imag(meij), choixIdx)
    
    frequence_elec = (omegaEtude*D)/(2*np.pi*bh[choixB]**2) # on revient à la fréquence thermique à partir de omega
    plt.semilogx(frequence_elec, deltaT_in_M)

plt.title("delta_T en fonction de la frequence thermique")
plt.xlabel("Frequence thermique")
plt.ylabel("Température moyenne (en °C)")
plt.show()

# V3omega
for choixB in range(0,1): # remplacer par range(bh.size) pour parcourir tous les b (attention c'est illisible)
    # extraction des valeurs qui vont nous intéresser
    choixIdx = df.index[(df["bh"] == bh[choixB])&(df["ts"]==ts[0])].to_numpy() # retourne les index des lignes qui vérifient b = bh[choixTs] et ts = ts[choixTs]
    omegaEtude = omega
    V3omega_in = np.real(V3omega)
    plt.loglog(omega, V3omega_in)


plt.title("V3omega en fonction de omega")
plt.xlabel("omega")
plt.ylabel("v3omega")
plt.show()

# amplitude et phase en fonction de la frequence
fig, ax1 = plt.subplots()
for choixB in range(0,1): # remplacer par range(bh.size) pour parcourir tous les b (attention c'est illisible)
    # extraction des valeurs qui vont nous intéresser
    choixIdx = df.index[(df["bh"] == bh[choixB])&(df["ts"]==ts[0])].to_numpy() # retourne les index des lignes qui vérifient b = bh[choixTs] et ts = ts[choixTs]
    omegaEtude = np.take(omega, choixIdx)
    ampEtude = np.take(amplitude, choixIdx)
    phaseEtude = np.take(phase, choixIdx)
    
    frequence_elec = (omegaEtude*D)/(2*np.pi*bh[choixB]**2) # on revient à la fréquence electrique à partir de omega
    ax1.semilogx(frequence_elec, ampEtude, color="red", label="amplitude")
    ax1.set_xlabel("frequence électrique")
    ax1.set_ylabel("amplitude")
    ax2 = ax1.twinx()
    ax2.semilogx(frequence_elec, phaseEtude, color="green", label="phase")
    ax2.set_ylabel("phase")

plt.title("amplitude et phase en fonction de la frequence electrique")
ax1.legend()
ax2.legend(loc="center right")
plt.show()



# retracer les asymptotes suivant le format defini ici ----> me faire un retour sur le tracé asymptotique
# ajouter la fonction MeijerG dans la methode f_u _ utiliser la librairie MeijerG comprendre les coefficients.
# integration numerique via Simpsom dans f_u ou separemment


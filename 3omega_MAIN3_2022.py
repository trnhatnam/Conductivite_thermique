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
frequence = np.arange(1, 5000, 200)          #domaine d`étude frequentiel

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
    V3omega_asympt = 0.5 * V0 * TCR * asympt    # calculate thrid harmonic from DT

    #ajouter fonctions exacxtes calculée en fonction de MeijerG ou Simpson

    #faire les calculs ici : amplitude , phase, V3omega

    return asympt, V3omega_asympt   #ajouter fle retour des calculs


# return a tuple of array. Remember to assign two otypes.
# f_u_vec = np.vectorize(f_u, otypes=[np.complex128,      # MeijerG is complex number ou Simpson
#                                     np.complex128,      # asympt is complex number
#                                     np.ndarray,         # amplitude is array
#                                     np.ndarray,         # phase is array
#                                     np.complex128,      # V3omega_asympt is complex number
#                                     np.complex128]      # V3omega is complex number issu de MeijerG ou Simpson
#                                     )

#tup = f_u_vec(omega) # tuple of arrays: (val1, asympt, amplitude, phase, V3omega_asympt, V3omega)

############## Fait par TRINH Nhat-nam
## partie calcul
f_u_vec2 = np.vectorize(f_u, otypes=[np.complex128, np.complex128]) # juste pour tester le traçage asymptotique
tup = f_u_vec2(omega) # maitenant, tup : (asympt, V3omega_asympt)

# extraction de la partie réel de l'asymptote
asympt = tup[0]
moduleT = np.abs(asympt)
argT = np.angle(asympt)
deltaT = moduleT*np.cos(argT)

# choix du tracé : on demande une valeur pour b et ts
print("Options pour b :", bh)
choixB = int(input("Avec quelle valeur de b voulez-vous travailler (tapez l'indice du b dans le tableau ci-dessus) ? : "))
print("Options pour ts :", ts)
choixTs = int(input("Avec quelle valeur de ts voulez-vous travailler (tapez l'indice de ts dans le tableau ci-dessus) ? : "))

# extraction des valeurs qui vont nous intéresser
choixIdx = df.index[(df["bh"] == bh[choixB])&(df["ts"]==ts[choixTs])].tolist() # retourne les index des lignes qui vérifient b = bh[choixTs] et ts = ts[choixTs]
omegaEtude = np.take(omega, choixIdx)
tempMoyEtude = np.take(deltaT, choixIdx)

## partie plot
def omg2lamb(omg):
    # omega -> Profondeur de pénétration (en µm)
    return np.sqrt((bh[choixB]**2)/(2*omg))/(1e-6)

fig, ax1 = plt.subplots()

# traçage de l'asymptote
plt.semilogx(omegaEtude, tempMoyEtude, linestyle="--", label="température moyenne") # traçage semilog en x
plt.title("Approximation de la température moyenne à basse fréquence en fonction de omega")
plt.xlabel("Omega")
plt.ylabel("Température moyenne (en °C)")

# ajout de la profondeur de pénétration comme deuxième abscisse x
ax2 = ax1.twiny()
ax2.invert_xaxis()
ax2.set_xlabel("Profondeur de pénétration en µm")
ax2.semilogx(omg2lamb(omegaEtude), tempMoyEtude)

ax1.legend()
plt.show()

# retracer les asymptotes suivant le format defini ici ----> me faire un retour sur le tracé asymptotique
# ajouter la fonction MeijerG dans la methode f_u _ utiliser la librairie MeijerG comprendre les coefficients.
# integration numerique via Simpsom dans f_u ou separemment
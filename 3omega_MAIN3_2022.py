from mpmath import meijerg, j
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import axes3d
%matplotlib qt # nécessaire pour Spyder, vous pouvez enlever sur Vscode


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

df['Thermal_freq'] = thermal_freq       # Insert thermal frequency into data frame
df['T_depth'] = T_depth                 # Insert thermal penetration depth into data frame
df["Re"] = np.real(tup[0])              # Insert val1 real part into data frame
df["Im"] = np.imag(tup[0])              # Insert val1 imagine part into data frame
df["Asympt_Re"] = np.real(tup[1])       # Insert asympt real part into data frame
df["Asympt_Im"] = np.imag(tup[1])       # Insert asympt imagine part into data frame
df['Amplitude'] = tup[2]                # Insert amplitude into data frame
df['Phase'] = tup[3]                    # Insert phase into data frame
df["V3asympt_Re"] = np.real(tup[4])     # Insert V3omega_asympt real part into data frame
df["V3asympt_Im"] = np.imag(tup[4])     # Insert V3omega_asympt imagine part into data frame
df["V3_Re"] = np.real(tup[5])           # Insert V3omega real part into data frame
df["V3_Im"] = np.imag(tup[5])           # Insert V3omega imagine part into data frame

# sauvegarde des resultats dans un fichier .csv generé selon bh
for bh_elem in bh:
    fname = f"bh={bh_elem:.4e}.csv"
    df_save = df[(df["bh"] == bh_elem)]
    df_save.to_csv(fname)

# température vs fréquence électrique vs t_depth
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('frequence électrique en Hz')
ax.set_ylabel('température en K')
ax.set_zlabel('t_depth en µm')
cl = ["black","blue","green"]
t = ts[0]
for b in bh:  
    # partie reel
    ax.plot(df[(df['bh'] == b)&(df['ts'] == t)]['frequence'],
            df[(df['bh'] == b)&(df['ts'] == t)]["Re"], 
            df[(df['bh'] == b)&(df['ts'] == t)]['T_depth'],
            linewidth=0.75, 
            zdir="z",
            label="b="+str(b),
            color = cl[np.where(bh == b)[0][0]])
    # partie imaginaire
    ax.plot(df[(df['bh'] == b)&(df['ts'] == t)]['frequence'],
                df[(df['bh'] == b)&(df['ts'] == t)]["Im"],
                df[(df['bh'] == b)&(df['ts'] == t)]['T_depth'],
                linewidth=0.75, 
                zdir="z",
                linestyle='--',
                color=cl[np.where(bh == b)[0][0]])

ax.set_title('Température vs fréquence électrique vs T_depth')
ax.legend()
plt.show()


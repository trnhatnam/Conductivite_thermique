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
 ###TEST SIMPSON

####

def simpson(f,a,b,N):
    ai = np.linspace(a,b,N+1)
    h = (b-a)/N
    somme=0
    for i in range (N):
        somme+=f(ai[i])+4*f((ai[i]+ai[i+1])/2)+f(ai[i+1])
    return (h/6)*somme

def f(y,freq):
    return ((np.sin(y)*np.sin(y))/y**2)*(1/(y**2 + j*freq))

P = V0**2/R0 #puissance par unité de longueur
bh = 1e-6
a = 10e-6
b = 10e6
N = 1000
#y_tab = np.linspace(a,b,N)
freq_tab= np.linspace(-10,10,N)


omega2 = (bh**2)*2* np.pi * frequence / D
Result = []
for freq in omega2 : #f pour fréquence

    Intregrale_freq = []
    y_tab = np.linspace(a,b,N+1)
    h = (b-a)/N
    simpson_result = 0
    for i in range (y_tab.size) :
         simpson_result+=f(y_tab[i],freq)+4*f((y_tab[i]+y_tab[i+1])/2,freq)+f(y_tab[i+1],freq)
    T_freq = P*(h/6)*simpson_result
    Result.append(T_freq)


def f_u(omega_elem):
     
    asympt = (P / (k*L*math.pi)) * (-(1 / 2) * np.log(omega_elem) + 3 / 2 - gamma - j * ((math.pi) / 4))

    #ajouter fonctions exacxtes calculée en fonction de MeijerG ou Simpson
    # meijerg
    val1 = (-j*P / (4*L*k*math.pi * omega_elem)) * meijerg([[1, 3 / 2], []], [[1, 1], [0.5, 0]], j * omega_elem)  #solution approximée via fnction MeijerG on recupere reel et imaginaire
    
    #simpson
    N = 100 #voir la meilleure valeur à prendre avec le prof
    valsimpson = (P/math.pi*k)*simpson(f,10*(-10),10*10,N) #N sous intervalles à définir
    #faire les calculs ici : amplitude , phase, V3omega
    amplitude = math.sqrt(np.real(val1) ** 2 + np.imag(val1) ** 2)
    phase = math.degrees(math.atan(np.imag(val1) / np.real(val1)))
    V3omega_asympt = 0.5 * V0 * TCR * asympt                                        # calculate thrid harmonic from DT
    V3omega = 0.5 * V0 * TCR * val1            

    return val1, valsimpson, asympt, amplitude, phase, V3omega, V3omega_asympt   #ajouter fle retour des calculs


#return a tuple of array. Remember to assign two otypes.
f_u_vec = np.vectorize(f_u, otypes=[np.complex128,      # MeijerG is complex number ou Simpson
                                    np.complex128,      # asympt is complex number
                                    np.ndarray,         # amplitude is array
                                    np.ndarray,         # phase is array
                                    np.complex128,      # V3omega_asympt is complex number
                                    np.complex128]      # V3omega is complex number issu de MeijerG ou Simpson
                                    )

tup = f_u_vec(omega) # tuple of arrays: (val1, asympt, amplitude, phase, V3omega_asympt, V3omega)
meij, asympt = tup[0], tup[1]

for choixB in range(0,1): # remplacer par range(bh.size) pour parcourir tous les b (attention c'est illisible)
    # extraction des valeurs qui vont nous intéresser
    choixIdx = df.index[(df["bh"] == bh[choixB])&(df["ts"]==ts[0])].to_numpy() # retourne les index des lignes qui vérifient b = bh[choixTs] et ts = ts[choixTs]
    omegaEtude = np.take(omega, choixIdx)
    deltaT_in = np.take(np.real(asympt), choixIdx)
    deltaT_out = np.take(np.imag(asympt), choixIdx)
    deltaT_in_M = np.take(np.real(meij), choixIdx)
    deltaT_out_M = np.take(np.imag(meij), choixIdx)

    # traçage de l'asymptote
    plt.semilogx(omegaEtude, deltaT_in, linestyle="--",linewidth=0.75,label="bf in phase b = " + str(bh[choixB])) # traçage semilog en x
    plt.semilogx(omegaEtude, deltaT_out, linestyle="--",linewidth=0.75, label="bf out phase b = " + str(bh[choixB]))
    
    # traçage de meijerg
    plt.semilogx(omegaEtude, deltaT_in_M, linewidth=1.25, label="meijerg in phase b = " + str(bh[choixB]))
    plt.semilogx(omegaEtude, deltaT_out_M, linewidth=1.25, label="meijerg out phase b = " + str(bh[choixB]))

plt.title("Asymptotes de la température moyenne à basse fréquence en fonction de omega")
plt.xlabel("Omega")
plt.ylabel("Température moyenne (en °C)")

plt.legend()
plt.show()
####
# retracer les asymptotes suivant le format defini ici ----> me faire un retour sur le tracé asymptotique
# ajouter la fonction MeijerG dans la methode f_u _ utiliser la librairie MeijerG comprendre les coefficients.
# integration numerique via Simpsom dans f_u ou separemment5

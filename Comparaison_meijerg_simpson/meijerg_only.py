#########################
# LIBRAIRIES A IMPORTER #
#########################
from mpmath import meijerg, j
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
from matplotlib import rc

# An instance of RcParams for handling default Matplotlib values.
# params = {
#     'figure.figsize': [16, 9],
#     'axes.labelsize' : 12,
#     'figure.autolayout': True,
#     'mathtext.rm': 'Times New Roman',
#     'font.family': 'Times New Roman',
#     'font.serif': "Times New Roman",
#     }
# matplotlib.rcParams.update(params)

# Resolution complex number
# mp.dps = 15                                 # Mpmath setting, dps is the decimal precision
# mp.pretty = True                            # Setting the mp.pretty option will use the str()-style output for repr() as well
# a = mpf(0.25)                               # a = 0.25
# b = mpf(0.25)                               # b = 0.25
# z = mpf(0.75)                               # z = 0.25

########################
# PARAMETRAGE DU MODELE#
########################

# PROPRIETES DU MATERIAU
D = 8.8e-5                               # Thermal diffusivity. 1.14e-3 diamond diffusivity, 200 -1000 e-8 for BN - CNT D = 4.6e-4, 0.025e4 polymers PVP aqueous, 1.3e-7
k = 140                                     # Thermal conductivity. Unit: W/mK (Diamond - 220 - 420 for BN) CNT 750 W/mK  or from 50 to 80 W/mK, 2000 W/mK, 0.27 PVP
gamma = 0.5772                              # Euler constant
V0 = 1                                      # Fundamental peak value
R0 = 80                                     # Room temp resistance of the heater
TCR = 0.00253                               # Temperature Coefficient of the electrical Resistivity (TRC), usually in the ~10^-3 /K. (Au TCR = 0.002535 /K)
L = 0.002                                   # Length of the heater
P = V0**2/R0                         # W/m power per unit length

# VARIABLES
start_f = -3
end_f = 7
points_f = 10000
frequency = np.logspace(start_f, end_f, points_f)                 # Frequency range f
# bh = np.array([5e-6, 10e-6, 15e-6, 20e-6])              # Heater half width bh
bh = np.array([5e-6])              # Heater half width bh

ts = np.array([350e-6])                                 # Substrate thickness ts
linear_limit = 25*D/(4*(math.pi)*(ts**2))               # Limits linear regimes
planar_limit = np.array(D/(100*(math.pi)*(bh**2)))      # Limits planar regimes

# # PRODUIT CARTESIEN DE VARIABLES POUR CREER UNE DATAFRAME df
idx = pd.MultiIndex.from_product([bh, ts,  frequency], names=["bh", 'ts', "frequency"])
df = pd.DataFrame(index=idx).reset_index()

omega = (df["bh"].values**2 * df["frequency"].values) * (2 * math.pi / D)           # Omega is vectorized naturally.
thermal_freq = (2 * df['frequency'])                                                # Thermal frequency is 2nd harmonic
T_depth = (np.sqrt(2*D / (2 * (math.pi) * thermal_freq)))/1e-6                      # Thermal penetration depth function of omega thermal not electrical - Cahill sqrt(2D/wt)

# METHODE DE RESOLUTION ANALYTIQUE PAR FONCTIONS MEIJERG
def f_u(omega_elem):
    # The function meijerg is evaluated as a combination of hypergeometric series
    val1 = (-j*P / (4*L*k*math.pi * omega_elem)) * meijerg([[1, 3 / 2], []], [[1, 1], [0.5, 0]], j * omega_elem)    # solution analytique fonction MeijerG on recupere reel et imaginaire
    asympt = (P / (k*L*math.pi)) * (-(1 / 2) * np.log(omega_elem) + 3 / 2 - gamma - j * ((math.pi) / 4))            # simplifcation asymptotique
    # extraction : amplitude , phase, V3omega
    amplitude = math.sqrt(np.real(val1) ** 2 + np.imag(val1) ** 2)
    phase = math.degrees(math.atan(np.imag(val1) / np.real(val1)))
    V3omega_asympt = 0.5 * V0 * TCR * asympt                                        # # issu du 3ieme harmonique
    V3omega = 0.5 * V0 * TCR * val1                                                 

    return val1, asympt, amplitude, phase, V3omega_asympt, V3omega

# # Creer un tuple qui contient l`ensemble des elements calculés
f_u_vec = np.vectorize(f_u, otypes=[np.complex128,      # val1 is complex number
                                    np.complex128,      # asympt is complex number
                                    np.ndarray,         # amplitude is array
                                    np.ndarray,         # phase is array
                                    np.complex128,      # V3omega_asympt is complex number
                                    np.complex128]      # V3omega is complex number
                                    )

tup = f_u_vec(omega)                    # tuple of arrays: (val1, asympt, amplitude, phase, V3omega_asympt, V3omega)

# stocker les elements du tuple dans une dataframe - cela permettra de tracer plus simplement
df['Thermal_freq'] = thermal_freq       # Insert thermal frequency into data frame
df["Re"] = np.real(tup[0])              # Insert val1 real part into data frame
df["Im"] = np.imag(tup[0])              # Insert val1 imagine part into data frame


# sauvegarde des resultats dans un fichier .csv generé (juste pour test)
df.to_csv("meijerg_only.csv")

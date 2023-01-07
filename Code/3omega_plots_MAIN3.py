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

########################
# PARAMETRAGE DU MODELE#
########################
# rechercher :  silicium (D = 8.8e-5, k=140 W/mk)
# PROPRIETES DU MATERIAU
D =  8.8e-5                         # Thermal diffusivity. 1.14e-3 diamond diffusivity, 200 -1000 e-8 for BN - CNT D = 4.6e-4, 0.025e4 polymers PVP aqueous, 1.3e-7
k = 140                             # Thermal conductivity. Unit: W/mK (Diamond - 220 - 420 for BN) CNT 750 W/mK  or from 50 to 80 W/mK, 2000 W/mK, 0.27 PVP, 
gamma = 0.5772                              # Euler constant
V0 = 1                                      # Fundamental peak value
R0 = 80                                     # Room temp resistance of the heater
TCR = 0.00253                               # Temperature Coefficient of the electrical Resistivity (TRC), usually in the ~10^-3 /K. (Au TCR = 0.002535 /K)
L = 0.002                                   # Length of the heater
P = V0**2/(R0)                            # W/m power per unit length

# VARIABLES
frequency = np.logspace(0, 6, num=12, endpoint=True)                 # Frequency range f
# bh = np.array([5e-6, 10e-6, 15e-6, 20e-6])              # Heater half width bh
bh = np.array([5e-6, 10e-6])              # Heater half width bh

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

# Code couleur pour les courbes
colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
columns = ['frequency', 'Thermal_freq']
xlabels = ['Frequency (Hz)', 'Thermal excitation frequency (Hz)']
linestyles = ['-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted']

# graduation mineures et majeures
locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1, ),numticks=100)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=np.arange(2, 10) * .1,numticks=100)
fig = plt.figure()
ax1 = fig.add_subplot(121)

########################
# Figure 1 #
########################

### ax1
ax1.tick_params(which="both", axis="both",direction="in")
ax1.minorticks_on()                                               # Display minor ticks on the axes.
ax1.set_ylabel('<T> (K)')                                         # Set the y-axis label
ax1.set_xlabel(xlabels[0])                                        # Set the bottom x-axis label
ax1.set_xscale('log')                                             # Set the bottom x-axis in log scale
ax1.xaxis.set_major_locator(locmaj)
ax1.xaxis.set_minor_locator(locmin)
ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax1.margins(x=0)
axs_top1 = ax1.twiny()                                            # Create a twin Axes sharing the y-axis.
axs_top1.invert_xaxis()                                                 # Invert the top x-axis.
axs_top1.tick_params(which="both", axis="both",direction="in")
axs_top1.set_xlabel('Thermal penetration depth: q$^\mathrm{-1}$ = (2D/\u03C9$_\mathrm{t}$)$^\mathrm{1/2}$' + ' (\u03BCm)', labelpad=10)
axs_top1.set_xticks(df['T_depth'])                                      # Set the top x-axis' tick locations.
axs_top1.set_xscale('log')                                              # Set the top x-axis in log scale
axs_top1.xaxis.set_major_locator(locmaj)
axs_top1.xaxis.set_minor_locator(locmin)
axs_top1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axs_top1.margins(x=0)

for axis in [ax1.xaxis, ax1.yaxis, axs_top1.xaxis, axs_top1.yaxis]:
    formatter = FuncFormatter(lambda x, _: '{:.5g}'.format(x))
    axis.set_major_formatter(formatter)


ax1.text(10, 0, r'Out-of-phase', fontsize=10)
ax1.text(10, np.mean(df[(df["bh"] == bh[0]) & (df["ts"] == ts[0])]['Re']), r'In-phase', fontsize=10)


for bh_array in bh:
    # Re and Imag part vs Frequency and Thermal penetration depth
    # Real part Vs Frequency, bh = to bh[1]
    column = columns[0]
    ax1.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])][column], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Re'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='bh=' + '{0:.1f}'.format((bh_array/1e-6)) + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{0:.0f}'.format((ts[0]/1e-6)) + 'um')
    # Real part Vs Thermal penetration depth, bh = to bh[1]
    axs_top1.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['T_depth'], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Re'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='bh=' + '{0:.1f}'.format((bh_array/1e-6))   + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{0:.0f}'.format((ts[0]/1e-6)) + 'um')

    # Imag part Vs Frequency, bh = to bh[1]
    ax1.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])][column], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Im'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='_nolegend_')
    # Imagine part Vs Thermal penetration depth, bh = to bh[1]
    axs_top1.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['T_depth'], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Im'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='_nolegend_')

# Asympt Real part Vs frequency, bh = to bh[1], choose the larger number from the bh array
ax1.plot(
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])][column], 
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['Asympt_Re'], 
    color= 'black', 
    linestyle='--',
    label='bh=' + '{0:.1f}'.format((bh[bh.size-1]/1e-6)) + '\u03BCm -  ' + r't$_{\rm s}$' + ' = ' + '{:.0f}'.format((ts[0]/1e-6)) + 'um  (asympt)')
# Asympt Imagine part Vs Thermal penetration depth, bh = to bh[1], choose the larger number from the bh array
ax1.plot(
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])][column], 
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['Asympt_Im'], 
    color= 'black', 
    linestyle='--',
    label='_nolegend_')
ax1.legend(loc='best', frameon=False, fontsize='8')


### ax2 :
    
ax2 = fig.add_subplot(122)
ax2.tick_params(which="both", axis="both",direction="in")
ax2.minorticks_on()                                               # Display minor ticks on the axes.
ax2.set_ylabel('<T> (K)')                                         # Set the y-axis label
ax2.set_xlabel(xlabels[1])                                        # Set the bottom x-axis label
ax2.set_xscale('log')                                             # Set the bottom x-axis in log scale
ax2.xaxis.set_major_locator(locmaj)
ax2.xaxis.set_minor_locator(locmin)
ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax2.margins(x=0)
axs_top2 = ax2.twiny()                                            # Create a twin Axes sharing the y-axis.
axs_top2.invert_xaxis()                                                 # Invert the top x-axis.
axs_top2.tick_params(which="both", axis="both",direction="in")
axs_top2.set_xlabel('Thermal penetration depth: q$^\mathrm{-1}$ = (2D/\u03C9$_\mathrm{t}$)$^\mathrm{1/2}$' + ' (\u03BCm)', labelpad=10)
axs_top2.set_xticks(df['T_depth'])                                      # Set the top x-axis' tick locations.
axs_top2.set_xscale('log')                                              # Set the top x-axis in log scale
axs_top2.xaxis.set_major_locator(locmaj)
axs_top2.xaxis.set_minor_locator(locmin)
axs_top2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axs_top2.margins(x=0)

for axis in [ax2.xaxis, ax2.yaxis, axs_top2.xaxis, axs_top2.yaxis]:
    formatter = FuncFormatter(lambda x, _: '{:.5g}'.format(x))
    axis.set_major_formatter(formatter)


ax2.text(10, 0, r'Out-of-phase', fontsize=10)
ax2.text(10, np.mean(df[(df["bh"] == bh[0]) & (df["ts"] == ts[0])]['Re']), r'In-phase', fontsize=10)


for bh_array in bh:
    # Re and Imag part vs Frequency and Thermal penetration depth
    # Real part Vs Frequency, bh = to bh[1]
    column = columns[1]
    ax2.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])][column], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Re'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='bh=' + '{0:.1f}'.format((bh_array/1e-6)) + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{0:.0f}'.format((ts[0]/1e-6)) + 'um')
    # Real part Vs Thermal penetration depth, bh = to bh[1]
    axs_top2.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['T_depth'], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Re'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='bh=' + '{0:.1f}'.format((bh_array/1e-6))   + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{0:.0f}'.format((ts[0]/1e-6)) + 'um')

    # Imag part Vs Frequency, bh = to bh[1]
    ax2.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])][column], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Im'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='_nolegend_')
    # Imagine part Vs Thermal penetration depth, bh = to bh[1]
    axs_top2.plot(
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['T_depth'], 
        df[(df["bh"] == bh_array) & (df["ts"] == ts[0])]['Im'], 
        colours[np.where(bh == bh_array)[0][0]], 
        label='_nolegend_')

# Asympt Real part Vs frequency, bh = to bh[1], choose the larger number from the bh array
ax2.plot(
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])][column], 
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['Asympt_Re'], 
    color= 'black', 
    linestyle="--",
    label='bh=' + '{0:.1f}'.format((bh[bh.size-1]/1e-6)) + '\u03BCm -  ' + r't$_{\rm s}$' + ' = ' + '{:.0f}'.format((ts[0]/1e-6)) + 'um  (asympt)')
# Asympt Imagine part Vs Thermal penetration depth, bh = to bh[1], choose the larger number from the bh array
ax2.plot(
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])][column], 
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['Asympt_Im'], 
    color= 'black', 
    linestyle="--",
    label='_nolegend_')
ax2.legend(loc='best', frameon=False, fontsize='8')



########################
# Figure 2 #
########################

### ax21
fig2 = plt.figure()
ax21 = fig2.add_subplot(121)

ax21.tick_params(which="both", axis="both",direction="in")
axs_top21 = ax21.twiny()
axs_top21.invert_xaxis()
axs_top21.tick_params(which="both", axis="both",direction="in")
ax21.minorticks_on()

# Plot the second figure
ax21.set_ylabel('Third harmonic voltage V3$_\mathrm{\u03C9}$ (V)')
ax21.set_xlabel('Thermal excitation frequency (Hz)')
ax21.set_xscale('log')
axs_top21.set_xticks(df['T_depth'])
axs_top21.set_xscale('log')
axs_top21.set_xlabel('Thermal penetration depth: q$^\mathrm{-1}$ = (2D/\u03C9$_\mathrm{t}$)$^\mathrm{1/2}$' + ' (\u03BCm)', labelpad=10)
ax21.xaxis.set_major_locator(locmaj)
axs_top21.xaxis.set_major_locator(locmaj)
ax21.xaxis.set_minor_locator(locmin)
axs_top21.xaxis.set_minor_locator(locmin)
ax21.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axs_top21.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
for axis in [ax21.xaxis, ax21.yaxis, axs_top21.xaxis, axs_top21.yaxis]:
    formatter = FuncFormatter(lambda x, _: '{:.5g}'.format(x))
    axis.set_major_formatter(formatter)
ax21.margins(x=0)
axs_top21.margins(x=0)
ax21.text(20, 0, r'Out-of-phase', fontsize=10)
ax21.text(20, np.mean(df[(df["bh"] == bh[0]) & (df["ts"] == ts[0])]['V3_Re']), r'In-phase', fontsize=10)

for bh_array2 in bh:
    # Re and Imag
    # V3 omega Real part Vs Thermal Frequency, bh = to bh[1]
    ax21.plot(
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['Thermal_freq'], 
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['V3_Re'], 
        colours[np.where(bh == bh_array2)[0][0]], 
        label='bh=' + '{:.1f}'.format((bh_array2/1e-6)) + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{:.0f}'.format((ts[0]/1e-6)) + 'um')
    # V3 omega Real part Vs Thermal penetration depth, bh = to bh[1]
    axs_top21.plot(
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['T_depth'], 
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['V3_Re'], 
        colours[np.where(bh == bh_array2)[0][0]], 
        label='bh=' + '{:.1f}'.format((bh_array2/1e-6)) + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{:.0f}'.format((ts[0]/1e-6)) + 'um')
    # V3 omega Imagine part Vs Thermal Frequency, bh = to bh[1]
    ax21.plot(
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['Thermal_freq'], 
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['V3_Im'], 
        colours[np.where(bh == bh_array2)[0][0]], 
        label='_nolegend_')
    # V3 omega Imagine part Vs Thermal penetration depth, bh = to bh[1]
    axs_top21.plot(
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['T_depth'], 
        df[(df["bh"] == bh_array2) & (df["ts"] == ts[0])]['V3_Im'], 
        colours[np.where(bh == bh_array2)[0][0]], 
        label='_nolegend_')
# V3asympt Real part Vs Thermal Frequency, bh = to bh[1]
ax21.plot(
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['Thermal_freq'], 
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['V3asympt_Re'], 
    color= 'black', 
    linestyle='--',  
    label='bh=' + '{0:.1f}'.format((bh[bh.size-1]/1e-6)) + '\u03BCm -  ' + r't$_{\rm s}$' + ' = ' + '{:.0f}'.format((ts[0]/1e-6)) + 'um (asympt)')
# V3asympt Imagine part Vs Thermal Frequency, bh = to bh[1]
ax21.plot(
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['Thermal_freq'], 
    df[(df["bh"] == bh[bh.size-1]) & (df["ts"] == ts[0])]['V3asympt_Im'], 
    color= 'black', 
    linestyle='--', 
    label='_nolegend_')

ax21.legend(loc='best', frameon=False, fontsize='8')


#### ax22
ax22 = fig2.add_subplot(122)
ax22.tick_params(which="both", axis="both",direction="in")
axs_top3 = ax22.twinx()
axs_top3.tick_params(which="both", axis="both",direction="in")
# The fourth figure setting
ax22.set_ylabel('Amplitude', fontsize=12)
ax22.yaxis.label.set_color('blue')
ax22.tick_params(axis='y', colors='blue')
axs_top3.set_ylabel('Phase (${^o}$)', fontsize=12)
axs_top3.yaxis.label.set_color('red')
axs_top3.tick_params(axis='y', colors='red')
ax22.set_xlabel('Frequency (Hz)')
ax22.set_xscale('log')
ax22.xaxis.set_major_locator(locmaj)
ax22.xaxis.set_minor_locator(locmin)
for axis in [ax22.xaxis, ax22.yaxis]:
    formatter = FuncFormatter(lambda x, _: '{:.16g}'.format(x))
    axis.set_major_formatter(formatter)
ax22.margins(x=0)
axs_top3.margins(x=0)

ax22.text(60, np.min(df[(df["bh"] == bh[0]) & (df["ts"] == ts[0])]['Amplitude'])-2, r'Amplitude', color="blue", fontsize=10)
ax22.text(60, np.mean(df[(df["bh"] == bh[0]) & (df["ts"] == ts[0])]['Phase']), r'Phase', color="red", fontsize=10)

for bh_array3 in bh:
    # Amplitude and phase
    # Amplitude Vs Frequency, bh = to bh[1]
    ax22.plot(
        df[(df["bh"] == bh_array3) & (df["ts"] == ts[0])]['frequency'], 
        df[(df["bh"] == bh_array3) & (df["ts"] == ts[0])]['Amplitude'], 
        color= 'blue', 
        linestyle=linestyles[np.where(bh == bh_array3)[0][0]],     
        label='bh=' + '{:.0f}'.format((bh_array3/1e-6)) + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{:.0f}'.format((ts[0]/1e-6)) + 'um')
    #Phase Vs Frequency, bh = to bh[1]
    ax22.plot(
        df[(df["bh"] == bh_array3) & (df["ts"] == ts[0])]['frequency'], 
        df[(df["bh"] == bh_array3) & (df["ts"] == ts[0])]['Phase'], 
        color= 'red', 
        linestyle=linestyles[np.where(bh == bh_array3)[0][0]],      
        label='bh=' + '{:.0f}'.format((bh_array3/1e-6)) + '\u03BCm - ' + r't$_{\rm s}$' + ' = ' + '{:.0f}'.format((ts[0]/1e-6)) + 'um')

ax22.legend(loc='best', frameon=False, fontsize='8')
plt.show()

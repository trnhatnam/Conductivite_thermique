import numpy as np
import cmath
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import matplotlib.ticker
import pandas as pd


#Integration de l`integrale complexe par la methode de Simspon
#le modele considere une integration sur [0, inf] mais il s`agit de considerer l`interval d`integration infini
# i.e. prendre un pas petit

# interval d`integration de l`integrale et definition du pas pour Simpson
y_0 = 1e-10
y_inf = 1   # suffisant à `mais prendre n grand
n = 10000 * y_inf
h = (y_inf - y_0) / (n-1) 
print("h parameter: ", h)
y_points = np.linspace(y_0, y_inf, n)  
print("Number of points: ", len(y_points))

#gamme de frequence
# start_f = 0.001
# end_f = 1e6
# step_f = 10000
# frequency = np.arange(start_f, end_f+step_f, step_f)   # 
start_f = -3
end_f = 7
points_f = 10000
frequency = np.logspace(start_f, end_f, points_f)
print("Number of frequency: ", len(frequency))
angular_frequency = 2 * (np.pi) * frequency
b = 5e-6
D = 8.8e-5
omega = ((b**2) * angular_frequency) / D  # puslation depend de la geometrie bh

#proprietes du matieriau
k = 140
V0 = 1
R0 = 80
P = V0**2/(R0)
L = 0.002
PD = P/(np.pi * k) # Power density

def f(freq,y):
    return (PD/L) * (np.sin(y)**2) / (y**2 * np.emath.sqrt(y**2 + freq * 1j))

#integration
simpson_points = np.zeros(omega.size)
for i in range(omega.size):
    freq = omega[i]
    print("Omega sweep: {:.4f}rad/s / {:.4f}rad/s, {:.4f}%".format(freq, max(omega), freq/max(omega)*100))
    simpson_points[i] = np.sum((h/6)*(f(freq,y_points[:-1]) + 4*(f(freq,(y_points[1:] + y_points[:-1])/2)) + f(freq,y_points[1:])))

# x = frequency
x = omega
y_real = np.real(simpson_points)
y_imag = np.imag(simpson_points)

df = pd.DataFrame()
df['Thermal_freq'] = frequency
df["Re"] = y_real
df["Im"] = y_imag
df.to_csv("simpson.csv")

plt.semilogx(frequency, y_real)
plt.show()


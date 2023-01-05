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
points_f = 100
frequency = np.logspace(start_f, end_f, points_f)
print("Number of frequency: ", len(frequency))
angular_frequency = 2 * (np.pi) * frequency
b = 5e-6
D = 2.67e-7
omega = ((b**2) * angular_frequency) / D  # puslation depend de la geometrie bh

#proprietes du matieriau
k = 0.55
V0 = 1
R0 = 80
P = V0**2/(R0)
PD = P/(np.pi * k) # Power density

#integration
simpon_points = []
for freq in omega:
    print("Omega sweep: {:.4f}rad/s / {:.4f}rad/s, {:.4f}%".format(freq, max(omega), freq/max(omega)*100))
    numeric_points = []
    for y in y_points:
        if y == 0:
            print("Zero in the denominator")
        else:
            T = PD * (np.sin(y)**2) / (y**2 * cmath.sqrt(y**2 + freq * 1j))
            # T = (np.sin(y)**2) / (y**2 * cmath.sqrt(y**2 + freq * 1j))
            numeric_points.append(T)
    simpson = (h/3) * (numeric_points[0] + 2*sum(numeric_points[2:n-2:2]) + 4*sum(numeric_points[1:n-1:2]) + numeric_points[n-1])
    simpon_points.append(simpson)

# x = frequency
x = omega
y_real = np.real(simpon_points)
y_imag = np.imag(simpon_points)

plt.semilogx(frequency, y_real)
plt.show()

df = pd.DataFrame()
df['Thermal_freq'] = frequency
df["Re"] = y_real
df["Im"] = y_imag
df.to_csv("simpson.csv")
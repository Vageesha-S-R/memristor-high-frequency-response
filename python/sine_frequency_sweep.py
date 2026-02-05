"""
Frequency sweep of TaOx memristor under sinusoidal input.
Reproduces high-frequency averaging behavior from:
Messaris – High Frequency Response of Non-Volatile Memristors
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ================= Device parameters =================
A, B = 1e-10, 1e-4
sigma_on, sigma_off = 0.45, 0.013
x_on, x_off = 0.06, 0.4
sigma_p, beta = 4e-5, 500
Gm, b, a = 0.025, 4.7, 7.2e-6

k = 50
V0 = 0.5


# ================= Model =================
def G(x, v):
    return Gm*x + a*np.exp(b*np.sqrt(abs(v))) * (1-x)

def g_plus(x, v):
    return B*np.sinh(v/sigma_on)*np.exp(-(x**2)/(x_on**2))*np.exp((G(x,v)*v**2)/sigma_p)

def g_minus(x, v):
    x_safe = max(x, 1e-6)
    return A*np.sinh(v/sigma_off)*np.exp(-(x_off**2)/(x_safe**2))*np.exp(-(1 + beta*G(x,v)*v**2))

def fk(v):
    return 1/(1+np.exp(-k*v))


# ================= Simulation =================
def simulate_frequency(f, zoom=False):

    def vin(t):
        return V0*np.sin(2*np.pi*f*t)

    def ode(t, x):
        v = vin(t)
        return [g_plus(x[0], v)*fk(v) + g_minus(x[0], v)*fk(-v)]

    T = 1/f
    dt = T/150

    t = np.arange(0, 0.2e-3, dt)

    sol = solve_ivp(
        ode, (0, 0.2e-3), [0.5],
        method="LSODA", t_eval=t,
        max_step=dt, rtol=1e-11, atol=1e-14
    )

    x = sol.y[0]

    mask = t > (t[-1]-T)
    mean_x = np.trapz(x[mask], t[mask]) / T
    print(f"{f/1e3:.1f} kHz → mean = {mean_x:.6f}")

    if zoom:
        mask_zoom = t > (t[-1]-3*T)
        plt.figure()
        plt.plot(t[mask_zoom]*1e3, x[mask_zoom], 'r')
        plt.axhline(mean_x, ls='--', color='k')
        plt.xlabel("Time (ms)")
        plt.ylabel("x(t)")
        plt.ticklabel_format(style='plain', useOffset=False)
        plt.title("100 MHz – last 3 periods")
        plt.grid()
        plt.show()

    return t, x


# ================= Main =================
freqs = [30e3, 150e3, 500e3, 100e6]
colors = ['tab:blue','tab:orange','tab:green','tab:red']

plt.figure()

for f, c in zip(freqs, colors):
    t, x = simulate_frequency(f, zoom=(f==100e6))
    label = f"{f/1e3:.0f} kHz" if f < 1e6 else "100 MHz"
    plt.plot(t*1e3, x, c, label=label)

plt.xlabel("Time (ms)")
plt.ylabel("State x(t)")
plt.legend()
plt.grid()
plt.title("State evolution for multiple frequencies")
plt.show()

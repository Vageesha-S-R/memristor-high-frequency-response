"""
Fading memory test:
Different initial states converge to the same steady state
under 100 MHz sinusoidal excitation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ================= Parameters =================
A, B = 1e-10, 1e-4
sigma_on, sigma_off = 0.45, 0.013
x_on, x_off = 0.06, 0.4
sigma_p, beta = 4e-5, 500
Gm, b, a = 0.025, 4.7, 7.2e-6

V0 = 0.5
f = 100e6
k = 50


# ================= Model =================
def vin(t):
    return V0*np.sin(2*np.pi*f*t)

def G(x, v):
    return Gm*x + a*np.exp(b*np.sqrt(abs(v))) * (1-x)

def g_plus(x, v):
    return B*np.sinh(v/sigma_on)*np.exp(-(x**2)/(x_on**2))*np.exp((G(x,v)*v**2)/sigma_p)

def g_minus(x, v):
    x_safe = max(x, 1e-6)
    return A*np.sinh(v/sigma_off)*np.exp(-(x_off**2)/(x_safe**2))*np.exp(-(1 + beta*G(x,v)*v**2))

def fk(v):
    return 1/(1+np.exp(-k*v))

def ode(t, x):
    v = vin(t)
    return [g_plus(x[0], v)*fk(v) + g_minus(x[0], v)*fk(-v)]


# ================= Simulation =================
x0_list = [0.2, 0.4, 0.6, 0.8, 1.0]

plt.figure()

for x0 in x0_list:

    sol = solve_ivp(
        ode, (0, 0.2e-3), [x0],
        method="LSODA",
        t_eval=np.linspace(0, 0.2e-3, 10000)
    )

    plt.plot(sol.t*1e3, sol.y[0], label=f"x0={x0}")

plt.xlabel("Time (ms)")
plt.ylabel("State x(t)")
plt.title("Fading memory (100 MHz sine)")
plt.legend()
plt.grid()
plt.show()

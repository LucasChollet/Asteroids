import numpy as np
import matplotlib.pyplot as plt

G = 2.95824e-4

# Sun
Ms = 1  # 1.988e30
mu_s = G * Ms
UA = 149.5978707

# Jupiter
Mj = 0.000950000000475  # 1 / 1047.348625
aJ = 5.2
Tj = 2 * np.pi * np.sqrt(np.power(aJ, 3) / G)
wJ = 2 * np.pi / Tj

m = 0
a = 8.25
test = [a, 0.0, 0.0, 0.0, np.sqrt(G) / np.sqrt(a), 0.0]


def rk4(y, t, step, f):
    half_step = step / 2
    k1 = f(y, t)
    k2 = f(y + half_step * k1, t + half_step)
    k3 = f(y + half_step * k2, t + half_step)
    k4 = f(y + step * k3, t + step)
    return y + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def integrate(y, t):
    return np.hstack((y[3:6], compute_gravity_3body(y[0:3], t)))


def compute_gravity_3body(v, t):
    uJ = np.array([aJ * np.cos(wJ * t), aJ * np.sin(wJ * t), 0])
    # base_sun = - G * Ms * m / np.power(np.linalg.norm(v), 3)
    # base_jup = - G * Mj * m / np.power(np.linalg.norm(v - uJ), 3)
    # return base_sun * v + base_jup * (v - uJ)
    return - G * (Ms + m) / np.power(np.linalg.norm(v), 3) * v - G * Mj * ((v - uJ) / np.power(np.linalg.norm(v - uJ), 3) + uJ / np.power(np.linalg.norm(uJ), 3))


def compute_semi_major_axis(v):
    return np.power((2 / np.linalg.norm(v[0:3]) - np.power(np.linalg.norm(v[3:6]), 2) / mu_s), -1)


def compute_eccentricity(v):
    return np.linalg.norm(np.cross(v[3:6], np.cross(v[0:3], v[3:6])) / mu_s - v[0:3] / np.linalg.norm(v[0:3]))


traj = plt.figure(1).add_subplot(111)
sma = plt.figure(2).add_subplot(111)
ecc = plt.figure(3).add_subplot(111)
traj.scatter(0, 0, c='r')
traj.axis("equal")

number_it = 365 * 30
step = 1
v = np.array(test)
print(v)
traj.scatter(v[0], v[1], c='C0', s=0.1)
for i in range(number_it):
    v = rk4(v, i * step, step, integrate)
    if i % 20 == 0:
        traj.scatter(v[0], v[1], c='C0', s=0.1)
        traj.scatter(aJ * np.cos(wJ * i * step), aJ * np.sin(wJ * i * step), c='C1', s=0.3)
        a = compute_semi_major_axis(v)
        e = compute_eccentricity(v)
        sma.scatter(i * step / 365, a, c='C1', s=0.1)
        ecc.scatter(i * step / 365, e, c='C2', s=0.1)

plt.show()
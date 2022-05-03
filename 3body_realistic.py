import numpy as np
import matplotlib.pyplot as plt

# Orbital element

def kep2cart(a, e, i, w, W, t, M0, t0):
    n = k / np.power(a, 3/2)
    M = M0 + n * (t - t0)
    # Newton
    E = M
    while E - e * np.sin(E) - M > 0.01:
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

    # On ne calcule pas l'anomalie vraie ?

    r = a * (1 - e * np.cos(E))

    X = a * (np.cos(E) - e)
    Y = a * np.sqrt(1 - np.power(e, 2)) * np.sin(E)
    X_d = - (n * np.power(a, 2) * np.sin(E)) / r
    Y_d = (n * np.power(a, 2) * np.sqrt(1 - np.power(e, 2)) * np.cos(E)) / r


    # Define rotation matrix
    R3W = create_rot3(-W)
    R1i = create_rot1(-i)
    R3w = create_rot3(-w)

    return R3W @ R1i @ R3w @ np.array([[X, Y, 0]]).T, R3W @ R1i @ R3w @ np.array([[X_d, Y_d, 0]]).T


def create_rot1(t):
    return np.array([[1, 0, 0], [0, np.cos(t), np.sin(t)], [0, -np.sin(t), np.cos(t)]])


def create_rot3(t):
    return np.array([[np.cos(t), np.sin(t), 0], [-np.sin(t), np.cos(t), 0], [0, 0, 1]])


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
    uJ = kep2cart(aJ, eJ, iJ, wJ, WJ, t, M0J, T0J)[:, 0]
    return - G * (Ms + mA) / np.power(np.linalg.norm(v), 3) * v - G * mJ * (
            (v - uJ) / np.power(np.linalg.norm(v - uJ), 3) + uJ / np.power(np.linalg.norm(uJ), 3))


def compute_semi_major_axis(v):
    return np.power((2 / np.linalg.norm(v[0:3]) - np.power(np.linalg.norm(v[3:6]), 2) / mu_s), -1)


def compute_eccentricity(v):
    return np.linalg.norm(np.cross(v[3:6], np.cross(v[0:3], v[3:6])) / mu_s - v[0:3] / np.linalg.norm(v[0:3]))


def compute_inclination(v):
    return np.arccos(((np.cross(v[0:3], v[3:6])) / (np.linalg.norm(v[0:3]) * np.linalg.norm(v[3:6])))[2])


G = 2.95824e-4
k = np.sqrt(G)

### Sun
Ms = 1  # 1.988e30
mu_s = G * Ms
UA = 149.5978707

### Jupiter
mJ = 1 / 1047.348625
aJ = 5.202575
eJ = 0.048908
iJ = 1.3038 * np.pi / 180
WJ = 100.5145 * np.pi / 180
wJ = 273.8752 * np.pi / 180
M0J = 80.0392 * np.pi / 180
T0J = 2456600.5
# Tj = 2 * np.pi * np.sqrt(np.power(aJ, 3) / G)
#
# wtJ = 2 * np.pi / Tj

### Asteroid 2007 VW266
mA = 0
aA = 5.454
eA = 0.3896
iA = 108.358 * np.pi / 180
wA = 226.107 * np.pi / 180
WA = 276.509
M0A = 146.88
T0A = T0J


traj = plt.figure(1).add_subplot(111)
sma = plt.figure(2).add_subplot(111)
ecc = plt.figure(3).add_subplot(111)
inc = plt.figure(4).add_subplot(111)
traj.scatter(0, 0, c='r')
traj.axis("equal")

number_it = 365 * 100
step = 1
v = kep2cart(aA, eA, iA, wA, WA, 0, M0A, T0A)[:, 0]
print(v)
traj.scatter(v[0], v[1], c='C0', s=0.1)
for i in range(number_it):
    v = rk4(v, i * step, step, integrate)
    if i % 20 == 0:
        traj.scatter(v[0], v[1], c='C0', s=0.1)
        uJ = kep2cart(aJ, eJ, iJ, wJ, WJ, i * step, M0J, T0J)
        traj.scatter(uJ[0], uJ[1], c='C1', s=0.3)
        a = compute_semi_major_axis(v)
        e = compute_eccentricity(v)
        incl = compute_inclination(v) * 180 / np.pi
        sma.scatter(i * step / 365, a, c='C1', s=0.1)
        ecc.scatter(i * step / 365, e, c='C2', s=0.1)
        inc.scatter(i * step / 365, incl, c='C3', s=0.1)

plt.show()
import numpy as np
import matplotlib.pyplot as plt

G = 2.95824e-4

Ms = 1  # 1.988e30
mu_s = G * Ms
UA = 149.5978707

m = 1
test = [1.0, 0.0, 0.0, 0.0, np.sqrt(G), 0.0]


def rk4(y, t, step, f):
    half_step = step / 2
    k1 = f(y, t)
    k2 = f(y + half_step * k1, t + half_step)
    k3 = f(y + half_step * k2, t + half_step)
    k4 = f(y + step * k3, t + step)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6


def integrate(y, t):
    return np.hstack((y[3:6], compute_gravity_2body(y[0:3])))


def compute_gravity_2body(v):
    base = - G * Ms * m / np.power(np.linalg.norm(v), 3)
    return base * v


def compute_semi_major_axis(v):
    return np.power((2 / np.linalg.norm(v[0:3]) - np.power(np.linalg.norm(v[3:6]), 2) / mu_s), -1)


def compute_eccentricity(v):
    return np.linalg.norm(np.cross(v[3:6], np.cross(v[0:3], v[3:6])) / mu_s - v[0:3] / np.linalg.norm(v[0:3]))


traj = plt.figure(1).add_subplot(111)
sma = plt.figure(2).add_subplot(111)
ecc = plt.figure(3).add_subplot(111)
traj.scatter(0, 0, c='r')
traj.axis("equal")

number_it = 3650
step = 0.1
v = np.array(test)
traj.scatter(v[0], v[1], c='C0', s=0.1)
for i in range(number_it):
    v = rk4(v, i * step, step, integrate)
    if i % 10 == 0:
        traj.scatter(v[0], v[1], c='C0', s=0.1)
        a = compute_semi_major_axis(v)
        e = compute_eccentricity(v)
        sma.scatter(i * step, a, c='C1', s=0.1)
        ecc.scatter(i * step, e, c='C2', s=0.1)
#
plt.show()

#
# v[3:6] = -v[3:6]
# for i in range(number_it):
#     v = rk4(v, i*step, step, integrate)

# v_init = np.array(test)
#
# d_final = np.linalg.norm(v - v_init)
#
# print(f"Différence entre les deux positions: {d_final * UA} millions km")
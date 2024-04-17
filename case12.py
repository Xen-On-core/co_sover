import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

import numpy as np

B = 0
C = 0


def Kef(phi):
    return phi ** 3


def Phi(x, z, t):
    return 0.1 + 0.1 * ((z - C * t) ** 2 / ((z - C * t) ** 2 + 1) + (x - B * t) ** 2 / ((x - B * t) ** 2 + 1))


N = 101
M = 101
K = 101

pa = 1
ps0 = 1
rhof = 2
rhos = 2600
tau0 = 1 / 31104000
k0 = 10 ** (-13)
alpha = 101368
L = 5000
H = 50
g = 9.81
rhof_g = rhof * g * H / alpha

tau = 1 / (K - 1)
hx = 1 / (N - 1)
hz = 1 / (M - 1)

alph = tau0 * L ** 2 / (k0 * alpha)
phi0 = [[0] * N for _ in range(M)]
for j in range(M):
    for i in range(N):
        phi0[j][i] = Phi(i * hx, j * hz, 0)

v = [0.5 * np.exp(-pow((20.0 * i * hx - 10.0), 6.0)) for i in range(N)]
pf = [[[0] * N for _ in range(M)] for _ in range(K)]
ps = [[[0] * N for _ in range(M)] for _ in range(K)]

integrate_H = [0 for _ in range(N)]
integrate_Z = [[0] * N for _ in range(M)]

for i in range(N):
    integrate_H[i] = sum([1 / Kef(phi0[j][i]) for j in range(M)])
    for j in range(M):
        integrate_Z[j][i] = sum([1 / Kef(phi0[jj][i]) for jj in range(j)])

for k in range(K):
    vt = np.exp(-(10 * tau * k - 0.5) ** 10)
    for j in range(M):
        for i in range(N):
            pf[k][j][i] = pa + rhof_g * (1 - j * hz) - alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])
            ps[k][j][i] = ps0 + alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])


X = np.arange(0, 1 + 1 / (N - 1), 1 / (N - 1))
Y = np.arange(0, 1 + 1 / (M - 1), 1 / (M - 1))
X, Y = np.meshgrid(X, Y)
Z = np.array(pf)
Zs = np.array(ps)

font = {'size': 14}
matplotlib.rc('font', **font)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
ax.set_xlabel("$x$")
ax.set_ylabel("$z$")

ax.set(title="$P_f$")
sf = ax.plot_surface(X, Y, Z[10], cmap=cm.viridis)
plt.savefig("./case12_img/pf1.png")
sf.remove()

sf = ax.plot_surface(X, Y, Z[K-1], cmap=cm.viridis)
plt.savefig("./case12_img/pfT.png")
sf.remove()

ax.set(title="$P_s$")
sf = ax.plot_surface(X, Y, Zs[10], cmap=cm.viridis)
plt.savefig("./case12_img/ps1.png")
sf.remove()

sf = ax.plot_surface(X, Y, Zs[K-1], cmap=cm.viridis)
plt.savefig("./case12_img/psT.png")


Zv1 = Z[10, 0:101, 50]

v = [0.4 * np.exp(-pow((20.0 * i * hx - 10.0), 6.0)) for i in range(N)]
for k in range(K):
    for j in range(M):
        for i in range(N):
            vt = np.exp(-(10 * tau * k - 0.5) ** 10)
            pf[k][j][i] = pa + rhof * g * (1 - j * hz) - alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])
            ps[k][j][i] = ps0 + alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])

Z = np.array(pf)
Zv2 = Z[10, 0:101, 50]

v = [0.3 * np.exp(-pow((20.0 * i * hx - 10.0), 6.0)) for i in range(N)]
for k in range(K):
    for j in range(M):
        for i in range(N):
            vt = np.exp(-(10 * tau * k - 0.5) ** 10)
            pf[k][j][i] = pa + rhof * g * (1 - j * hz) - alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])
            ps[k][j][i] = ps0 + alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])

Z = np.array(pf)
Zv3 = Z[10, 0:101, 50]

v = [0.2 * np.exp(-pow((20.0 * i * hx - 10.0), 6.0)) for i in range(N)]
for k in range(K):
    for j in range(M):
        for i in range(N):
            vt = np.exp(-(10 * tau * k - 0.5) ** 10)
            pf[k][j][i] = pa + rhof * g * (1 - j * hz) - alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])
            ps[k][j][i] = ps0 + alph * phi0[j][i] * v[i] * vt * (integrate_Z[j][i] - integrate_H[i])

Z = np.array(pf)
Zv4 = Z[10, 0:101, 50]

fig, ax = plt.subplots()
sf.remove()

fig, ax = plt.subplots()
X1 = np.arange(0, 1 + 1 / (N - 1), 1 / (N - 1))
ax.set_xlabel("$z$")
ax.set_ylabel("$P_f$")
ax.plot(X1, Zv1, color='blue', label="$v=0.5$")
ax.plot(X1[::4], Zv1[::4], marker='v')
ax.plot(X1, Zv2, color='red', label="$v=0.4$")
ax.plot(X1[::4], Zv2[::4], marker='*')
ax.plot(X1, Zv3, color='green', label="$v=0.3$")
ax.plot(X1[::4], Zv3[::4], marker='o')
ax.plot(X1, Zv4, color='orange', label="$v=0.2$")
ax.plot(X1[::4], Zv4[::4], marker='+')
plt.legend()
ax.grid()
plt.savefig("./case12_img/pf1_v.png")

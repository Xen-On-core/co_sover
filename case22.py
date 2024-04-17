import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

import numpy as np


def B(x, t):
    return 0.5 #* 15 * (x-0.5) / (1 - pow(abs(15 * (x-0.5)), 2))


def C(x, t):
    return 0.5


def Kef(phi):
    return phi ** 3


N = 51
M = 51
K = 60001

pa = 1
ps0 = 1
rhof = 2
rhos = 2600
tau0 = 1 / 691200
k0 = 10 ** (-13)
alpha = 101368
L = 2000
H = 50
g = 9.81
rhof_g = rhof * g * H / alpha

tau = tau0 #1 / (K - 1)
hx = 1 / (N - 1)
hz = 1 / (M - 1)

alph = tau0 * L ** 2 / (k0 * alpha)
phi0 = [[[0.82] * N for _ in range(M)] for _ in range(K)]
phi0_12 = [[[0.82] * N for _ in range(M)] for _ in range(K)]

alphx = [0 for _ in range(N)]
betax = [0 for _ in range(N)]
alphz = [0 for _ in range(M)]
betaz = [0 for _ in range(M)]
for k in range(K - 1):
    for j in range(1, M-1):
        alphx[0] = 0
        betax[0] = 0.82
        for i in range(N-1):
            BB = B((i + 1) * hx, tau * k) - B((i - 1) * hx, tau * k)
            A = G = B(i * hx, tau * k) * tau / (4 * hx)
            P = 1
            F = phi0[k][j][i] - C(i * hx, tau * k) * tau / (4 * hz) * (phi0[k][j+1][i] - phi0[k][j-1][i]) - phi0[k][j][i] * tau / (4*hx) * BB
            alphx[i+1] = G / (P - alphx[i]*A)
            betax[i+1] = (F + A*betax[i]) / (P - alphx[i]*A)

        phi0_12[k+1][j][N-1] = 0.82 #betax[N-1]/(1 - alphx[N-1])
        for i in range(N-2, -1, -1):
            phi0_12[k + 1][j][i] = alphx[i+1] * phi0_12[k+1][j][i+1] + betax[i+1]

    for i in range(1, N-1):
        alphz[0] = 0
        betaz[0] = 0.82
        for j in range(M-1):
            BB = B((i + 1) * hx, tau * k) - B((i - 1) * hx, tau * k)
            A = G = C(i * hx, tau * k) * tau / (4 * hz)
            P = 1
            F = phi0_12[k][j][i] - B(i * hx, tau * k) * tau / (4 * hx) * (phi0_12[k][j][i+1] - phi0_12[k][j][i-1]) - phi0[k][j][i] * tau / (4*hx) * BB
            alphz[j+1] = G / (P - alphz[j] * A)
            betaz[j+1] = (F + A * betaz[j]) / (P - alphz[j] * A)

        phi0[k+1][M-1][i] = betaz[M-1]/(1 - alphz[M-1])
        for j in range(M-2, -1, -1):
            phi0[k + 1][j][i] = alphz[j+1] * phi0[k+1][j+1][i] + betaz[j+1]

v = [0.5 * np.exp(-pow((20.0 * i * hx - 10.0), 10.0)) for i in range(N)]
pf = [[[0] * N for _ in range(M)] for _ in range(K)]
ps = [[[0] * N for _ in range(M)] for _ in range(K)]

integrate_H = [[0] * N for _ in range(K)]
integrate_H2 = [[0] * N for _ in range(K)]
integrate_Z = [[[0] * N for _ in range(M)] for _ in range(K)]
integrate_Z2 = [[[0] * N for _ in range(M)] for _ in range(K)]

ps_integration_H = [[0] * N for _ in range(K)]
ps_integration_Z = [[[0] * N for _ in range(M)] for _ in range(K)]

for k in range(K):
    for i in range(N):
        integrate_H[k][i] = sum([hz / Kef(phi0[k][j][i]) for j in range(M)])
        integrate_H2[k][i] = sum([j * hz / Kef(phi0[k][j][i]) for j in range(M)])
        for j in range(M):
            integrate_Z[k][j][i] = sum([hz / Kef(phi0[k][jj][i]) for jj in range(j)])
            integrate_Z2[k][j][i] = sum([jj * hz / Kef(phi0[k][jj][i]) for jj in range(j)])

for k in range(K):
    vt = np.exp(-(10 * tau * k - 0.5) ** 10)
    for j in range(M):
        for i in range(N):
            BB = (B((i + 1) * hx, tau * k) - B((i - 1) * hx, tau * k)) / (2 * hx)
            pf[k][j][i] = pa + rhof * g * (1 - j * hz) + alph * BB * (integrate_Z2[k][j][i] - integrate_H2[k][i]) - \
                          alph * (phi0[k][j][i] * v[i] * vt + BB) * (integrate_Z[k][j][i] - integrate_H[k][i])

X = np.arange(0, 1 + hx, hx)
Y = np.arange(0, 1 + hz, hz)
X, Y = np.meshgrid(X, Y)
Z = np.array(pf)
# Zs = np.array(ps)
Zphi = np.array(phi0)

font = {'size': 10}
matplotlib.rc('font', **font)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
ax.set_xlabel("$x$")
ax.set_ylabel("$z$")

ax.set(title="$\phi$")
sf = ax.plot_surface(X, Y, Zphi[10], cmap=cm.viridis)
plt.savefig("./case22_img/phi1.png")
sf.remove()

sf = ax.plot_surface(X, Y, Zphi[K - 1], cmap=cm.viridis)
plt.savefig("./case22_img/phiT.png")
sf.remove()

ax.set(title="$P_f$")
sf = ax.plot_surface(X, Y, Z[10], cmap=cm.viridis)
plt.savefig("./case22_img/pf1.png")
sf.remove()

sf = ax.plot_surface(X, Y, Z[K - 1], cmap=cm.viridis)
plt.savefig("./case22_img/pfT.png")
sf.remove()

# ax.set(title="$P_s$")
# sf = ax.plot_surface(X, Y, Zs[10], cmap=cm.viridis)
# plt.savefig("./case22_img/ps1.png")
# sf.remove()
#
# sf = ax.plot_surface(X, Y, Zs[K - 1], cmap=cm.viridis)
# plt.savefig("./case22_img/psT.png")
# sf.remove()

# ax.set(title="$P_{tot}$")
# sf = ax.plot_surface(X, Y, Zphi[10] * Z[10] - (1 - Zphi[10]) * Zs[10], cmap=cm.viridis)
# plt.savefig("./case22_img/ptot1.png")
# sf.remove()
#
# sf = ax.plot_surface(X, Y, Zphi[K - 1] * Z[K - 1] - (1 - Zphi[K - 1]) * Zs[K - 1], cmap=cm.viridis)
# plt.savefig("./case22_img/ptotT.png")

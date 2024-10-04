import numpy as np
from matplotlib import pyplot as plt

# T = np.logspace(-6, 2, 120)
# BC = np.exp(1 / T)
# BC_approx1 = 1 + 1 / T
# eps = 1e-2
# BC_approx2 = 1 + 1 / np.sqrt(T * T + eps * eps)
# plt.plot(T, BC, color="black", label="BC")
# plt.loglog()
# plt.plot(T, BC_approx1, color="red", label="BC approx")
# plt.plot(T, BC_approx2, color="blue", label="BC approx eps e-2")

# #plt.vlines(x,())
# plt.ylim((1e-1, 1e10))
# plt.legend()
# plt.savefig("bcs_approx.png")


# w = laplacian(phi)
# need exp(3-phi/T)~1
# phi = 3T


def sigmoid(t, k=0.1, ft=0.5):
    tmax = np.max(t) * ft
    ts = (t - tmax / 2) / tmax
    return 1 / (1 + np.exp(-ts / k))


nsteps = 10000
tmax = 100.0
t = np.linspace(0.0, tmax, nsteps)

for k, ft in [(0.07, 0.5), (0.07, 0.75), (0.07, 0.9)]:
    plt.plot(t, sigmoid(t, ft=ft, k=k), label=f"k = {k}, ft={ft}")
plt.legend()
plt.savefig("sigmoid.png")

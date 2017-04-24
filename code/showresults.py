import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def logsumexp(xs):
    """
    logsumexp.
    """
    biggest = xs.max()
    return np.log(np.sum(np.exp(xs - biggest))) + biggest

# Load results, normalise prior weights
sample_info = pd.read_csv("sample_info.csv")
logps = sample_info["logX"] - logsumexp(sample_info["logX"])

def canonical(temperatures):
    """
    Construct a canonical distribution.
    """
    assert len(temperatures) == sample_info.columns.size - 2

    temp = logps.copy()
    for i in range(0, len(temperatures)):
        temp += sample_info.iloc[:,i+2] / temperatures[i]
    logZ = logsumexp(temp)
    post = np.exp(temp - logZ)
    H = np.sum(post * (temp - logZ - logps))

    return {"logZ": logZ, "H": H}


x = np.linspace(0.0, 1.0, 20001)
def truth(T1, T2):
    p = np.exp(-(x - 0.5)**2/T1 - np.abs(x)/T2)
    Z = np.trapz(p, x=x)
    H = np.trapz(p/Z*np.log(p/Z + 1E-300), x=x)

    logZ = 10*np.log(Z)
    H *= 10

    return {"logZ": logZ, "H": H}

# Calculate log(Z) and H for some canonical distributions
T1 = 10.0**(np.linspace(-2.0, 4.0, 101))
T2 = T1.copy()
[T1, T2] = np.meshgrid(T1, T2[::-1])
logZ = T1.copy()
logZ_est = T1.copy()
H = T1.copy()
H_est = T1.copy()

for i in range(0, T1.shape[0]):
    for j in range(0, T1.shape[1]):
        temp1 = truth(T1[i, j], T2[i, j])
        temp2 = canonical([T1[i, j], T2[i, j]])
        logZ[i, j] = temp1["logZ"]
        logZ_est[i, j] = temp2["logZ"]
        H[i, j] = temp1["H"]
        H_est[i, j] = temp2["H"]
    print(i+1, "/", T1.shape[0])

# Plot points.
plt.plot(sample_info["scalars[0]"], sample_info["scalars[1]"],
         "k.", markersize=1, alpha=0.2)
plt.xlabel(r"$L_1$")
plt.ylabel(r"$L_2$")
plt.show()

# Plot phase diagrams
plt.figure(figsize=(9, 7))
plt.subplot(2, 2, 1)
extent = np.log10(np.array([T1.min(), T1.max(), T2.min(), T2.max()]))
plt.imshow(logZ_est, extent=extent)
plt.ylabel(r"$\log_{10} T_2$")
plt.title(r"$\log(Z)$")

plt.subplot(2, 2, 2)
plt.imshow(H_est, extent=extent)
plt.title(r"$H$")

plt.subplot(2, 2, 3)
error = logZ_est - logZ
plt.imshow(error, cmap="coolwarm", extent=extent,
           vmin=-np.max(np.abs(error)), vmax=np.max(np.abs(error)))
plt.xlabel(r"$\log_{10} T_1$")
plt.ylabel(r"$\log_{10} T_2$")

plt.subplot(2, 2, 4)
error = H_est - H
plt.imshow(error, cmap="coolwarm", extent=extent,
           vmin=-np.max(np.abs(error)), vmax=np.max(np.abs(error)))
plt.xlabel(r"$\log_{10} T_1$")
plt.show()


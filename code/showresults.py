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

# Plot the logL-logX curves of the individual scalars, then
# remove them
for i in [0, 1]:
    try:
        which = sample_info["which_scalar"] == i
        logX = sample_info["logX"][which]
        s = np.array(sample_info["scalars[{i}]".format(i=i)][which])
        ymin = s[len(s) // 20]
        ymax = s[-1]
        yrange = ymax - ymin

        plt.subplot(2, 1, i+1)
        plt.plot(logX, s, "-")
        plt.ylim([ymin, ymax + 0.05*yrange])
        plt.xlabel(r"$\log X$")
        plt.ylabel(r"$S_{k}$".format(k=i))
    except:
        pass
plt.show()


run0 = sample_info[sample_info["which_scalar"] == 0]
run1 = sample_info[sample_info["which_scalar"] == 1]
run2 = sample_info[sample_info["which_scalar"] == 2]
logps = run2["logX"] - logsumexp(run2["logX"])
depth = -(run2["logX"].min())



# Plot points.
plt.plot(run0["scalars[0]"], run0["scalars[1]"],
         "b.", markersize=1, alpha=0.1)
plt.plot(run1["scalars[0]"], run1["scalars[1]"],
         "r.", markersize=1, alpha=0.1)
plt.plot(run2["scalars[0]"], run2["scalars[1]"],
         "k.", markersize=1, alpha=0.1)
plt.xlabel(r"$S_1$")
plt.ylabel(r"$S_2$")
plt.show()

def canonical(temperatures):
    """
    Construct a canonical distribution.
    """
    assert len(temperatures) == run2.columns.size - 3

    temp = logps.copy()
    for i in range(0, len(temperatures)):
        temp += run2.iloc[:,i+3] / temperatures[i]
    logZ = logsumexp(temp)
    post = np.exp(temp - logZ)
    H = np.sum(post*(temp - logZ - logps))

    S1 = np.sum(post*run2["scalars[0]"])
    S2 = np.sum(post*run2["scalars[1]"])

    return {"logZ": logZ, "H": H, "S1": S1, "S2": S2}


x = np.linspace(0.0, 1.0, 20001)
def truth(T1, T2):
    p = np.exp(-(x - 0.5)**2/T1 - np.abs(x)/T2)
    Z = np.trapz(p, x=x)
    H = np.trapz(p/Z*np.log(p/Z + 1E-300), x=x)

    logZ = 10*np.log(Z)
    H *= 10

    return {"logZ": logZ, "H": H}

# Calculate log(Z) and H for some canonical distributions
T1 = 10.0**(np.linspace(-2.0, 2.0, 51))
T2 = T1.copy()
[T1, T2] = np.meshgrid(T1, T2[::-1])
logZ = T1.copy()
logZ_est = T1.copy()
H = T1.copy()
H_est = T1.copy()
S1_expected = T1.copy()
S2_expected = T2.copy()

for i in range(0, T1.shape[0]):
    for j in range(0, T1.shape[1]):
        temp1 = truth(T1[i, j], T2[i, j])
        temp2 = canonical([T1[i, j], T2[i, j]])
        logZ[i, j] = temp1["logZ"]
        logZ_est[i, j] = temp2["logZ"]
        H[i, j] = temp1["H"]
        H_est[i, j] = temp2["H"]
    print(i+1, "/", T1.shape[0])


# Plot phase diagrams
plt.figure(1, figsize=(9, 7))

cmap = plt.cm.viridis
cmap.set_bad("white", 1.0)
bad = H_est > 0.8*depth

plt.subplot(2, 2, 1)
extent = np.log10(np.array([T1.min(), T1.max(), T2.min(), T2.max()]))

plt.imshow(np.ma.array(logZ_est, mask=bad), extent=extent, cmap=cmap)
plt.ylabel(r"$\log_{10} T_2$")
plt.title(r"$\log(Z)$")

plt.subplot(2, 2, 2)
plt.imshow(np.ma.array(H_est, mask=bad), extent=extent, cmap=cmap)
plt.title(r"$H$")

plt.subplot(2, 2, 3)

plt.imshow(np.ma.array(S1_expected, mask=bad),
           extent=extent,
           cmap=cmap)
plt.xlabel(r"$\log_{10} T_1$")
plt.ylabel(r"$\log_{10} T_2$")
plt.title(r"$\left<S_1\right>$")

plt.subplot(2, 2, 4)
plt.imshow(np.ma.array(S2_expected, mask=bad),
           extent=extent,
           cmap=cmap)
plt.ylabel(r"$\log_{10} T_2$")
plt.title(r"$\left<S_2\right>$")

plt.figure(2)

cmap = plt.cm.coolwarm
cmap.set_bad("white", 1.0)

plt.subplot(1, 2, 1)
error = np.ma.array(logZ_est - logZ, mask=bad)
plt.imshow(error, cmap=cmap, extent=extent,
           vmin=-np.max(np.abs(error)), vmax=np.max(np.abs(error)))
plt.xlabel(r"$\log_{10} T_1$")
plt.ylabel(r"$\log_{10} T_2$")

plt.subplot(1, 2, 2)
error = np.ma.array(H_est - H, mask=bad)
plt.imshow(error, cmap=cmap, extent=extent,
           vmin=-np.max(np.abs(error)), vmax=np.max(np.abs(error)))
plt.xlabel(r"$\log_{10} T_1$")
plt.show()


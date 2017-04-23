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

    return {"logZ": logZ}

# Calculate normalising constant of a canonical distribution
print("logZ(T1=0.1, T2=1.0) = {lz}"\
            .format(lz=canonical([0.1, 1.0])["logZ"]))

# Plot points.
plt.plot(sample_info["scalars[0]"], sample_info["scalars[1]"],
         "b.", markersize=1, alpha=0.2)
plt.show()


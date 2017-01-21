import numpy as np
import matplotlib.pyplot as plt

output = np.loadtxt("sample_info.txt")

plt.plot(output[:,1], output[:,2], "k.", alpha=0.2)
plt.show()


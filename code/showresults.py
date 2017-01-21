import numpy as np
import matplotlib.pyplot as plt

output = np.loadtxt("output.txt")

plt.scatter(output[:,0], output[:,1], color="g", marker="o", s=output[:,2],
                alpha=0.2)
plt.show()


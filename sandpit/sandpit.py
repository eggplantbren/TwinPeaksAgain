import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

# Initial particles
N = 1
x = rng.rand(N)
y = 1.0 - x # rng.rand(N)

def calculate_ucc(xx, yy):
    ucc = 0
    for i in range(0, N):
        if x[i] > xx and y[i] > yy:
            ucc += 1
    return ucc

x_grid = np.linspace(0.0, 1.0, 1001)
y_grid = np.linspace(0.0, 1.0, 1001)
[x_grid, y_grid] = np.meshgrid(x_grid, y_grid)
y_grid = y_grid[::-1, :]

ucc_grid = np.zeros((1001, 1001), dtype="int64")
for i in range(0, 1001):
    for j in range(0, 1001):
        ucc_grid[i, j] = calculate_ucc(x_grid[i, j], y_grid[i, j])
    print(i+1)

# Plot the UCC map
plt.imshow(ucc_grid, interpolation="nearest",
                     extent=[0.0, 1.0, 0.0, 1.0], cmap="Blues")
plt.contour(x_grid, y_grid, ucc_grid, np.arange(0, N), colors="yellow")


# Plot the particles
plt.plot(x, y, "ro")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.show()


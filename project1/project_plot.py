import numpy as np
import matplotlib.pyplot as plt

#Read and plot all temperature grid
for t in range(10000):
 filename = f"temperature_{t}.txt"
 temperature_grid = np.loadtxt(filename)
 plt.imshow(temperature_grid, cmap='hot', interpolation='nearest')
 plt.colorbar()
 plt.title(f'Time Step {t}')
 plt.xlabel('X')
 plt.ylabel('Y')
 plt.show()

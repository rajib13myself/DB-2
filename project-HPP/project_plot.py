import numpy as np
import matplotlib.pyplot as plt

# Load temperature data from file
temperature_data = np.loadtxt("temperature_data.txt")

# Plot the temperature data
plt.imshow(temperature_data, cmap='hot', interpolation='nearest')
plt.colorbar(label='Temperature')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Temperature Distribution')
plt.show()

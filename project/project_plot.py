import numpy as np
import matplotlib.pyplot as plt

# Load temperature data from file
temperature_data = np.loadtxt("temperature_data.txt")

# Flip the temperature data upside down
temperature_data = np.flipud(temperature_data)

# Plot the temperature data
plt.imshow(temperature_data, cmap='hot', interpolation='nearest', origin='lower')
plt.colorbar(label='Temperature')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Temperature Distribution')
plt.show()

import matplotlib.pyplot as plt
import numpy as np

# Create a figure and axis
fig, ax = plt.subplots(figsize=(8, 1))

# Generate a range of values from 0 to 1
values = np.linspace(0, 1, 1000)

# Create a plasma colormap
cmap = plt.get_cmap('plasma')

# Plot the colormap
ax.imshow([values], cmap=cmap, aspect='auto')

# Remove axis ticks and labels
ax.set_xticks([])
ax.set_yticks([])

# Save the figure as a PNG
plt.savefig('plasma_colormap.png', bbox_inches='tight', dpi=150)

# Show the plot
plt.show()
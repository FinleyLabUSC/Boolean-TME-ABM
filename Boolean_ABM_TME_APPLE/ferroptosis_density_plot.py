import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ferroptosis_trigger_start = 15  # days until ferroptosis initiated
col_names = ['x_loc', 'y_loc', 'radius', 'state', 'pdl1', 'ferroptosis_sensitive']
BIN_WIDTH = 20 #um
TUMOR_EDGE = 1200 #um

# Calculate proportions for each set
proportions_all_sets = []

for set_idx in range(101):

    dir_path = f'dev_test/set_{set_idx}/cellLists/'  # location
    timestep_fpath = [f'{dir_path}day_{ferroptosis_trigger_start}_hour_{i}/cells.csv' for i in range(2, 25)]
    
    dfs = [pd.read_csv(fpath, names=col_names) for fpath in timestep_fpath[1:6]]

    for df in dfs:
        df['distance'] = np.sqrt(df['x_loc']**2 + df['y_loc']**2)
    
    bins = np.arange(0, TUMOR_EDGE, BIN_WIDTH)
    proportions = []

    for df in dfs:
        timestep_proportions = []
        for i in range(len(bins)-1):
            bin_mask = (df['distance'] >= bins[i]) & (df['distance'] < bins[i+1])
            bin_df = df[bin_mask]
            ferroptosis_sensitive_cells = bin_df[bin_df['ferroptosis_sensitive'] == -1]
            proportion = len(ferroptosis_sensitive_cells) / len(bin_df) if len(bin_df) > 0 else 0
            timestep_proportions.append(proportion)
        
        proportions.append(timestep_proportions)

    proportions_all_sets.append(proportions)

# Calculate mean proportions across all sets
mean_proportions = np.mean(proportions_all_sets, axis=0)

# Create 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('white')

for i, timestep_proportions in enumerate(mean_proportions):
    ax.plot(bins[:-1], timestep_proportions, zs=i, zdir='y', color=plt.cm.viridis(i/len(proportions)))
    

ax.set_zlim([0,1])
ax.set_xlim([0, TUMOR_EDGE])

ax.set_yticks(np.arange(0, 5))

ax.set_yticklabels(['1', '2', '3', '4', '5'], fontsize=10)

ax.set_xlabel('Distance from center (μm)', fontsize=10, fontweight='bold')
ax.set_ylabel('Timepoint', fontsize=10, fontweight='bold')
ax.set_zlabel('Proportion of \nferroptotic cells', fontsize=10, fontweight='bold')

ax.tick_params(axis='x', labelsize=8)
ax.tick_params(axis='y', labelsize=8)
ax.tick_params(axis='z', labelsize=8)
plt.savefig('ferroptosis_density_evolution', dpi=400)
plt.show()













import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np 
import random

rep_id = random.randint(0, 101)

day_number = 15 #same as ferroptosis initiation
hours = np.arange(2, 7)


fig, axs = plt.subplots(1, 5)

print(len(hours))
for i in range(len(hours)):
	

	fp = f'dev_test/set_{rep_id}/cellLists/day_{day_number}_hour_{hours[i]}/cells.csv'



	df = pd.read_csv(fp, header=None)
	df.columns = ["cellType", "x_loc", "y_loc","radius", "state", "pdl1_conc", "ferroptosis_sensitive"]
	df = df[df["cellType"] == 0] #only cancer cells

	sensitive_cells = df[df["ferroptosis_sensitive"] == -1]
	other_cells = df[df["ferroptosis_sensitive"] != -1]

	# Plot sensitive cells
	p = sns.scatterplot(data=sensitive_cells, ax = axs[i],x='x_loc', y='y_loc', color='cyan', alpha=0.9, s=1, legend=False)

	# Plot other cells
	sns.scatterplot(data=other_cells, ax = axs[i], x='x_loc', y='y_loc', hue='ferroptosis_sensitive', alpha=0.4, s=1, palette='plasma',legend=False)



	# Add title and labels
	xlim = p.get_xlim()
	ylim = p.get_ylim()
	p.set_xlim(-3500, 3500)
	p.set_ylim(-3500, 3500)
	p.set(xticklabels=[])
	p.set(yticklabels=[])
	p.set_facecolor('black')
	plt.xlabel('')
	plt.ylabel('')
	plt.tick_params(axis='x', which='both', bottom=False, top=False)
	plt.tick_params(axis='y', which='both', left=False, right=False)

plt.show()
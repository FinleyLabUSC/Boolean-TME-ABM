import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np 
import random


day_number = 12

fp = f'dev_test/set_0/cellLists/day_{day_number}/cells.csv'

fig = plt.figure(figsize=(6, 6))
df = pd.read_csv(fp, header=None)
df.columns = ["cellType", "x_loc", "y_loc","radius", "state", "pdl1_conc", "ferroptosis_sensitive"]
df = df[df["cellType"] == 0] #only cancer cells

p = sns.scatterplot(data=df, x='x_loc', y='y_loc', hue='ferroptosis_sensitive', alpha=0.9, s=1, legend=True)


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
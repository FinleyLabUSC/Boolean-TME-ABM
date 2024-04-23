import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation

# Function to update the plot for each day
def update_plot(day_number):
    plt.cla()  # Clear the previous plot
    
    fp = f'dev_test/set_0/cellLists/day_{day_number}/cells.csv'
    df = pd.read_csv(fp, header=None)
    df.columns = ["cellType", "x_loc", "y_loc","radius", "state", "pdl1_conc", "ferroptosis_sensitive"]
    df = df[df["cellType"] == 0] # Only cancer cells
    
    p = sns.scatterplot(data=df, x='x_loc', y='y_loc', hue='ferroptosis_sensitive', alpha=0.9, s=1, palette='plasma')
    p.get_legend().remove()


    
    # Add title and labels
    p.set_xlim(-3500, 3500)
    p.set_ylim(-3500, 3500)
    p.set(xticklabels=[])
    p.set(yticklabels=[])
    p.set_facecolor('black')
    # p.figure.colorbar(sm)
    plt.xlabel('')
    plt.ylabel('')
    plt.tick_params(axis='x', which='both', bottom=False, top=False)
    plt.tick_params(axis='y', which='both', left=False, right=False)
    
    plt.title(f'Day {day_number}')  # Set title
    
# Create the animation
fig = plt.figure(figsize=(6, 6))
ani = FuncAnimation(fig, update_plot, frames=range(25), interval=500)  # Iterate over day_number from 1 to 30 with a 500ms interval
ani.save('ferroptosis_evolution_test.mp4')
plt.show()
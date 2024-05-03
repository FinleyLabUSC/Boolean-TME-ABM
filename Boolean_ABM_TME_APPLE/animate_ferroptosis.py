import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches


# Function to update the plot for each day
def update_plot(day_number):
    plt.cla()  # Clear the previous plot
    
    fp = f'dev_test/set_0/cellLists/day_{day_number}/cells.csv'
    df = pd.read_csv(fp, header=None)
    df.columns = ["cellType", "x_loc", "y_loc","radius", "state", "pdl1_conc", "ferroptosis_sensitive"]
    df = df[df["cellType"] == 0] # Only cancer cells
    

    sensitive_cells = df[df["ferroptosis_sensitive"] == -1]
    other_cells = df[df["ferroptosis_sensitive"] != -1]

    
    # Plot sensitive cells
    ax = sns.scatterplot(data=sensitive_cells, x='x_loc', y='y_loc', color='cyan', alpha=0.9, s=1, legend=False)
    
    # Plot other cells
    sns.scatterplot(data=other_cells, x='x_loc', y='y_loc', hue='ferroptosis_sensitive', alpha=0.4, s=1, palette='plasma', ax=ax,legend=False)

    # p = sns.scatterplot(data=df, x='x_loc', y='y_loc', hue='ferroptosis_sensitive', alpha=0.9, s=1, palette='plasma')
    # plt.get_legend().remove()
    scale_length = 1000  # Length of the scale bar in micrometers
    scale_loc = (-3000, -3000)  # Location of the scale bar (bottom-left corner)
    scale_width = 100  # Width of the scale bar
    scale_height = 20  # Height of the scale bar
    scale_label_offset = 50  # Offset for the scale label
    plt.gca().add_patch(patches.Rectangle(scale_loc, scale_length, scale_height, color='white'))
    plt.text(scale_loc[0] + scale_length / 2, scale_loc[1] - scale_label_offset-200, f'{scale_length} μm', color='white', ha='center', fontsize=10)





    # Add title and labels
    plt.xlim(-3500, 3500)
    plt.ylim(-3500, 3500)
    plt.xticks([])
    plt.yticks([])
    plt.gca().set_facecolor('black')
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
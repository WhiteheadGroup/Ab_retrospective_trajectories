import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import os

plt.clf()

file1 = '/Users/siobhan/PycharmProjects/MAGMA-seq_proccessing_3/50thpctKds.csv'

file1_df = pd.read_csv(file1)

# Define unique categories
categories = ["002-S21F2", "1-57", "2-15", "C118", "CC12.1", "CC6.31", "4A8"]

# Define a colormap
# Define a list of 7 unique hexadecimal color codes
hex_colors = ['#FFC305', '#45ACB3', '#F27F34', '#3954B2', '#E24249', '#A451C2', '#9CDE62']

# Create a custom colormap
custom_cmap = ListedColormap(hex_colors)

# Plot each category with its assigned color from the custom colormap
plt.figure(figsize=(10, 6))

# Initialize a list for custom legend handles
legend_handles = []

for i, category in enumerate(categories):
    xdata = file1_df['Step']
    ydata = file1_df[category]

    # Plot scatter points
    plt.scatter(xdata, ydata, color=custom_cmap(i), marker='o', edgecolors=custom_cmap(i),
                facecolors='none', s=100, linewidths=2)

    # Plot lines

    # Plot lines with longer dashes
    line, = plt.plot(xdata, ydata, color=custom_cmap(i), linestyle='--', linewidth=0.75)
    line.set_dashes([10, 5])  # Set custom dash pattern

    # Create a custom legend handle combining both scatter and line
    custom_handle = Line2D([0], [0], color=custom_cmap(i), marker='o', markersize=8,
                           markerfacecolor='none', linestyle='--', linewidth=1.5)
    legend_handles.append((custom_handle, category))

# Unpack legend handles and labels
handles, labels = zip(*legend_handles)

plt.ylim(10, 1900)
plt.xlim(-0.075, 6.25)
plt.yscale('log')  # Set y-axis to logarithmic scale

# Set tick label sizes using Axes.tick_params()
ax = plt.gca()
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)

# Define font properties for axis labels and legend
font_properties = {'family': 'sans-serif', 'weight': 'normal', 'size': 14}

# Add labels, title, and legend with custom font
plt.xlabel('Development Step', fontdict=font_properties)
plt.ylabel('50th Percentile Kd >', fontdict=font_properties)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_position(('zero'))  # Place left spine at x=0

# Add custom legend
plt.legend(handles=handles, labels=labels, framealpha=1, prop=font_properties, frameon=False)

# Save the plot as a PNG file
#plt.savefig('50thPercentileKds.png', dpi=1200, bbox_inches='tight')

plt.show()
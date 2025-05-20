#!/usr/bin/env python3
# F3 Statistic 
# Figure 2D (varying alpha)
# This script computes an f3 statistic defined as (p3 – p1)·(p3 – p5) over a grid of p1 and p5 values, 
# then generates contour plots of this statistic for several fixed p3 values, arranging them in a 2×4 grid with a shared custom‐colored legend.

# Williams & Huber 2025 – Genomic Footprints 

import numpy as np                               # Numerical operations
import matplotlib.pyplot as plt                  # Plotting library
from matplotlib.colors import LinearSegmentedColormap  # For creating custom colormaps

def f3_statistic(p3, p1, p5):
    """
    Compute the f3 statistic:
      f3 = (p3 - p1) * (p3 - p5)
    This measures the product of allele frequency differences.
    """
    return (p3 - p1) * (p3 - p5)

def plot_f3_statistic(ax, p3):
    # 1. Generate p1 and p5 values on a 100×100 grid
    p1_values = np.linspace(0, 1, 100)    # p1 from 0 → 1
    p5_values = np.linspace(1, 0, 100)    # p5 reversed: 1 → 0 for visual orientation
    p1_grid, p5_grid = np.meshgrid(p1_values, p5_values)

    # 2. Evaluate the statistic over the grid
    f3_grid = f3_statistic(p3, p1_grid, p5_grid)

    # 3. Create a custom blue–white–red colormap
    colors = ['darkblue', 'blue', 'white', 'red', 'darkred']
    n_bins = 500
    cmap = LinearSegmentedColormap.from_list("custom", colors, N=n_bins)

    # 4. Set fixed color scale limits to center the white midpoint
    vmin, vmax = -0.25, 0.25

    # 5. Draw filled contours of f3_grid
    contourf = ax.contourf(
        p1_grid, p5_grid, f3_grid,
        levels=50, cmap=cmap, vmin=vmin, vmax=vmax
    )
    # 6. Overlay the zero contour in black
    ax.contour(p1_grid, p5_grid, f3_grid, levels=[0], colors='k', linewidths=1)

    # 7. Add a dotted line marking the p3 coordinate on the plot
    ax.plot([0, p3, p3, 1], [p3, p3, 1, 1], 'k:', linewidth=2)

    # 8. Label axes and title with larger fonts for readability
    ax.set_xlabel('p1', fontsize=20)
    ax.set_ylabel('p5', fontsize=20)
    ax.set_title(f'f3 statistic for p3 = {p3:.2f}', fontsize=20)

    # 9. Increase tick label size
    ax.tick_params(axis='both', which='major', labelsize=18)

    return contourf  # Return the contour set for legend creation

# List of p3 values to visualize
p3_values = [0.01, 0.05, 0.1, 0.2, 0.25, 0.35, 0.45, 0.5]

# Determine subplot layout: 2 rows, compute required columns
num_rows = 2
num_cols = int(np.ceil(len(p3_values) / num_rows))

# Create the main figure and subplots grid
fig, axs = plt.subplots(
    num_rows, num_cols,
    figsize=(5 * num_cols, 5 * num_rows)
)
fig.suptitle('f3(p3-p1)(p3-p5) for different p3 values', fontsize=24)

# Flatten the axes array for easy iteration
axs_flat = axs.flatten()

# Generate a contour plot for each p3 value
for ax, p3 in zip(axs_flat, p3_values):
    contourf = plot_f3_statistic(ax, p3)

# Remove any unused subplot panels
for ax in axs_flat[len(p3_values):]:
    ax.remove()

# Adjust layout to make space for the super-title
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Create a separate figure solely for the colorbar legend
fig_legend = plt.figure(figsize=(4, 8))
ax_legend = fig_legend.add_subplot(111)

# Add the shared colorbar to the legend figure
cbar = fig_legend.colorbar(contourf, cax=ax_legend, extend='both')
cbar.set_label('f3(p3-p1)(p3-p5)', fontsize=24, labelpad=20)

# Increase colorbar tick label size and force ticks on the left
cbar.ax.tick_params(axis='both', which='major', labelsize=20)
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.yaxis.set_tick_params(labelleft=True)

# Define and label specific tick positions for clarity
cbar.set_ticks([-0.25, -0.125, 0, 0.125, 0.25])
cbar.set_ticklabels(['-0.25', '-0.125', '0', '0.125', '0.25'])

# Display both the contour plots and the legend
plt.show()

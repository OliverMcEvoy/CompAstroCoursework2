import pandas as pd
import matplotlib.pyplot as plt
import time

# Load data
start_time = time.time()
rayleigh_blue_bins_df = pd.read_csv(
    "rayleigh_blue_position_bins.csv",
    names=["bin_start", "bin_end", "x Coordinate", "y Coordinate", "z Coordinate", "q"],
)
print(f"Time to load rayleigh_blue_bins_df: {time.time() - start_time:.2f} seconds")

start_time = time.time()
rayleigh_other_bins_df = pd.read_csv(
    "rayleigh_other_position_bins.csv",
    names=["bin_start", "bin_end", "x Coordinate", "y Coordinate", "z Coordinate", "q"],
)
print(f"Time to load rayleigh_other_bins_df: {time.time() - start_time:.2f} seconds")

# Create a list of DataFrames and their corresponding titles
binned_dataframes = [rayleigh_blue_bins_df, rayleigh_other_bins_df]
titles = ["Blue Light", "Other Colours"]

# Colors for each subplot
colors = [
    ["blue", "purple", "navy"],
    ["orange", "red", "brown"],
]  # Different for each dataset

# Set up the figure (3 rows, 2 columns)
plt.figure(figsize=(10, 9))

# Define the metrics to plot
metrics = ["x Coordinate", "y Coordinate", "z Coordinate"]
metric_labels = [
    "Normalised distribution",
    "Normalised Distribution",
    " Normalised Distribution",
]

# Plot each binned dataset
for idx, (df, title) in enumerate(zip(binned_dataframes, titles)):
    for j, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = plt.subplot(3, 2, j * 2 + idx + 1)  # 3 rows, 2 columns
        plt.scatter(
            df["bin_start"], df[metric], color=colors[idx][j], s=15, alpha=0.8
        )  # Only scatter, different colors
        plt.xlabel(f"{metric}")

        # Add ylabel only to the leftmost plots (idx == 0)
        if idx == 0:
            plt.ylabel(f"{label}")
        if j == 0:
            plt.title(f"{title}")

        # Set different x-axis limits for "avg_z" plots
        if j == 2:  # Third plot in each row
            plt.xlim(-200, 2000)
        else:
            if idx == 0:
                plt.xlim(-300, 300)
            else:
                plt.xlim(-10, 10)

plt.tight_layout()
plt.savefig("binned_position_scatter.png")
plt.show()

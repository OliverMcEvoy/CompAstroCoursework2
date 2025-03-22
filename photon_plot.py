import pandas as pd
import matplotlib.pyplot as plt
import time

# Load the new binned position data with timing
start_time = time.time()
photon_bins_df = pd.read_csv(
    "photon_position_bins.csv",
    names=["bin_start", "bin_end", "avg_x", "avg_y", "avg_z", "q"],
)
print(f"Time to load photon_bins_df: {time.time() - start_time:.2f} seconds")

start_time = time.time()
rayleigh_blue_bins_df = pd.read_csv(
    "rayleigh_blue_position_bins.csv",
    names=["bin_start", "bin_end", "avg_x", "avg_y", "avg_z", "q"],
)
print(f"Time to load rayleigh_blue_bins_df: {time.time() - start_time:.2f} seconds")

start_time = time.time()
rayleigh_other_bins_df = pd.read_csv(
    "rayleigh_other_position_bins.csv",
    names=["bin_start", "bin_end", "avg_x", "avg_y", "avg_z", "q"],
)
print(f"Time to load rayleigh_other_bins_df: {time.time() - start_time:.2f} seconds")

# Create a list of DataFrames and their corresponding titles
binned_dataframes = [photon_bins_df, rayleigh_blue_bins_df, rayleigh_other_bins_df]
titles = ["Photons", "Blue Light", "Other Colours"]

# Set up the figure for the binned data
plt.figure(figsize=(20, 15))  # Larger figure to accommodate 3 rows and 4 columns

# Define the metrics to plot
metrics = ["avg_x", "avg_y", "avg_z", "q"]
metric_labels = ["Average x", "Average y", "z", "Average q"]

# Plot each binned dataset
for idx, (df, title) in enumerate(zip(binned_dataframes, titles)):
    for j, (metric, label) in enumerate(zip(metrics, metric_labels)):
        plt.subplot(3, 4, idx * 4 + j + 1)  # 3 rows, 4 columns
        plt.bar(df["bin_start"], df[metric], width=1.0, color="blue", alpha=0.6)
        plt.xlabel("z (Position)")
        plt.ylabel(label)
        plt.xlim(-2000, 2000)  # Adjust based on your data range

plt.tight_layout()
plt.savefig("binned_position_analysis.png")
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import time

# Load data with timing
start_time = time.time()
photons_df = pd.read_csv("photon.csv", names=["y", "x", "theta", "iterations"])
print(f"Time to load photons_df: {time.time() - start_time:.2f} seconds")

start_time = time.time()
rayleigh_other_df = pd.read_csv(
    "rayleigh_other.csv", names=["y", "x", "theta", "iterations"]
)
print(f"Time to load rayleigh_other_df: {time.time() - start_time:.2f} seconds")

start_time = time.time()
rayleigh_blue_df = pd.read_csv(
    "rayleigh_blue.csv", names=["y", "x", "theta", "iterations"]
)
print(f"Time to load rayleigh_blue_df: {time.time() - start_time:.2f} seconds")

# Create a list of DataFrames and their corresponding titles
dataframes = [photons_df, rayleigh_blue_df, rayleigh_other_df]
titles = ["Photons", "Blue light", "Other Colours"]

# Set up the figure with subplots for each dataset
fig, axes = plt.subplots(3, 4, figsize=(20, 15))
fig.suptitle("Photon and Rayleigh Scattering Analysis", fontsize=16)

# Iterate over each DataFrame and create its subplots
for idx, (df, title) in enumerate(zip(dataframes, titles)):
    # Plot 1: Density function of theta
    axes[idx, 0].hist(df["theta"], bins=50, color="blue", alpha=0.6)
    axes[idx, 0].set_xlabel("Theta (radians)")
    axes[idx, 0].set_ylabel("Density")

    # Plot 2: Histogram of iterations
    axes[idx, 1].hist(df["iterations"], bins=50, color="green", alpha=0.6)
    axes[idx, 1].set_xlabel("Iterations")
    axes[idx, 1].set_ylabel("Count")

    # Plot 3: Histogram of y values
    axes[idx, 2].hist(df["y"], bins=50, color="orange", alpha=0.7)
    axes[idx, 2].set_xlabel("y")
    axes[idx, 2].set_ylabel("Count")
    axes[idx, 2].set_xlim(-1500, 1500)

    # Plot 4: Histogram of x values (New plot using X column)
    axes[idx, 3].hist(df["x"], bins=50, color="purple", alpha=0.7)
    axes[idx, 3].set_xlabel("x")
    axes[idx, 3].set_ylabel("Count")
    axes[idx, 3].set_xlim(-1500, 1500)  # Adjust based on your data range

plt.tight_layout()
plt.savefig("photon_analysis.png")
plt.show()

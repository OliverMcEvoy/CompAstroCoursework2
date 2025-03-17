import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("photon_scattering_results.csv", names=["y", "theta", "iterations"])

# Set up the figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Density function of theta
axes[0, 0].hist(df["theta"], bins=10, color="blue", alpha=0.6)
axes[0, 0].set_title("Density Function of Theta (Angle)")
axes[0, 0].set_xlabel("Theta (radians)")
axes[0, 0].set_ylabel("Density")

# Plot 2: Histogram of iterations per photon
axes[0, 1].hist(df["iterations"], bins=20, color="green", alpha=0.6)
axes[0, 1].set_title("Distribution of Iterations per Photon")
axes[0, 1].set_xlabel("Iterations")
axes[0, 1].set_ylabel("Count")

# Plot 3
axes[1, 0].hist(df["y"], bins=20, color="orange", alpha=0.6)
axes[1, 0].set_title("Distribution of y values per Photon")
axes[1, 0].set_xlabel("y")
axes[1, 0].set_ylabel("Count")

# Plot 4: Line plot of iterations over index
axes[1, 1].scatter(df["y"], df["theta"], color="purple", alpha=0.7)
axes[1, 1].set_title("Intensity vs y ")
axes[1, 1].set_xlabel("Index")
axes[1, 1].set_ylabel("Iterations")

# Adjust layout
plt.tight_layout()
plt.show()

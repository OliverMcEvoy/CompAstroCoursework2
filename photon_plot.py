import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("photon_scattering_results.csv", names=["y", "theta", "iterations"])

# Set up the figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Density function of theta
axes[0, 0].hist(df["theta"], bins=30, density=True, color="blue", alpha=0.6)
axes[0, 0].set_title("Density Function of Theta (Angle)")
axes[0, 0].set_xlabel("Theta (radians)")
axes[0, 0].set_ylabel("Density")

# Plot 2: Histogram of iterations per photon
axes[0, 1].hist(df["iterations"], bins=20, color="green", alpha=0.6)
axes[0, 1].set_title("Distribution of Iterations per Photon")
axes[0, 1].set_xlabel("Iterations")
axes[0, 1].set_ylabel("Count")

# Plot 3: Scatter plot of y vs. theta
axes[1, 0].scatter(df["theta"], df["y"], alpha=0.7, color="red")
axes[1, 0].set_title("Photon Displacement vs. Angle")
axes[1, 0].set_xlabel("Theta (radians)")
axes[1, 0].set_ylabel("Y Displacement")

# Plot 4: Line plot of iterations over index
axes[1, 1].scatter(df["theta"], df["iterations"], color="purple", alpha=0.7)
axes[1, 1].set_title("Iterations per Photon over Time")
axes[1, 1].set_xlabel("Index")
axes[1, 1].set_ylabel("Iterations")

# Adjust layout
plt.tight_layout()
plt.show()

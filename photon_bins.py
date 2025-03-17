import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV file
data = pd.read_csv(
    "binned_intensity.csv",
    header=None,
    names=["mu", "intensity", "bin_edge_lower", "bin_edge_upper"],
)

# Extract mu, intensity, and bin edges
mu = data["mu"]
intensity = data["intensity"]
bin_edge_lower = data["bin_edge_lower"]
bin_edge_upper = data["bin_edge_upper"]

# Plot settings
fig, axes = plt.subplots(1, 1, figsize=(10, 6))  # Set figure size

# Plot the intensity as a function of mu
axes.scatter(mu, intensity, marker="o", color="red", label="Intensity")

# Shade the regions of each bin
for i in range(len(bin_edge_lower)):
    axes.axvspan(
        np.cos(bin_edge_lower[i]),
        np.cos(bin_edge_upper[i]),
        color="gray",
        alpha=0.2,
        label=f"Bin {i + 1}" if i == 0 else "",
    )

# Add labels and title
axes.set_xlabel(r"$\mu = \cos\theta$", fontsize=14)
axes.set_ylabel("Normalized Intensity", fontsize=14)
plt.title("Intensity as a Function of $\mu = \cos\theta$", fontsize=16)


# Add legend
axes.legend()

# Save the plot as a PNG file
plt.savefig("intensity_plot_with_bins.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()

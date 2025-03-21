import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV files into separate DataFrames
photon_data = pd.read_csv(
    "photon_bins.csv",
    header=None,
    names=["mu", "intensity", "bin_edge_lower", "bin_edge_upper"],
)

rayleigh_blue_data = pd.read_csv(
    "rayleigh_blue_bins.csv",
    header=None,
    names=["mu", "intensity", "bin_edge_lower", "bin_edge_upper"],
)

rayleigh_other_data = pd.read_csv(
    "rayleigh_other_bins.csv",
    header=None,
    names=["mu", "intensity", "bin_edge_lower", "bin_edge_upper"],
)

# Create a list of DataFrames for easier iteration
dataframes = [photon_data, rayleigh_blue_data, rayleigh_other_data]
titles = ["Photon Bins", "Rayleigh Blue Bins", "Rayleigh Other Bins"]

# Plot settings
fig, axes = plt.subplots(1, 3, figsize=(18, 6))  # 1 row, 3 columns for subplots
fig.suptitle("Intensity as a Function of $\mu = \cos\theta$", fontsize=16)

# Iterate over each DataFrame and plot
for idx, (data, ax, title) in enumerate(zip(dataframes, axes, titles)):
    # Extract mu, intensity, and bin edges
    mu = data["mu"]
    intensity = data["intensity"]
    bin_edge_lower = data["bin_edge_lower"]
    bin_edge_upper = data["bin_edge_upper"]

    # Adjust mu values for indices >= 5
    for i in range(len(mu)):
        if i >= 5:
            mu[i] = 2 - mu[i]

    # Plot the intensity as a function of mu
    ax.scatter(mu, intensity, marker="o", color="red", label="Intensity")

    # Shade the regions of each bin
    for i in range(len(bin_edge_lower)):
        if i >= 5:
            ax.axvspan(
                2 - np.cos(bin_edge_lower[i]),
                2 - np.cos(bin_edge_upper[i]),
                color="gray",
                alpha=0.2,
                label=f"Bin {i + 1}" if i == 0 else "",
            )
        ax.axvspan(
            np.cos(bin_edge_lower[i]),
            np.cos(bin_edge_upper[i]),
            color="gray",
            alpha=0.2,
            label=f"Bin {i + 1}" if i == 0 else "",
        )

    # Add labels and title for each subplot
    ax.set_xlabel(r"$\mu = \cos\theta$", fontsize=12)
    ax.set_ylabel("Normalized Intensity", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend()

# Adjust layout and save the plot
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make room for the main title
plt.savefig("intensity_plot_with_bins_subplots.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()

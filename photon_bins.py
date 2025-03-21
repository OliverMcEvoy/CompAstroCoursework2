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

# Create a 2-row, 3-column subplot layout
fig, axes = plt.subplots(2, 3, figsize=(18, 12))  # 2 rows, 3 columns
fig.suptitle("Intensity as a Function of $\mu$ and $\\theta$", fontsize=18)

# Iterate over each DataFrame and plot
for idx, (data, title) in enumerate(zip(dataframes, titles)):
    # Extract mu, intensity, and bin edges
    mu = data["mu"]
    intensity = data["intensity"]
    bin_edge_lower = data["bin_edge_lower"]
    bin_edge_upper = data["bin_edge_upper"]

    # Compute theta in degrees
    theta = np.arccos(mu)

    # Top row: Intensity vs mu
    ax_mu = axes[0, idx]
    ax_mu.scatter(mu, intensity, marker="o", color="red", label="Intensity")

    # Shade the bin regions
    for i in range(len(bin_edge_lower)):
        ax_mu.axvspan(
            bin_edge_lower[i],
            bin_edge_upper[i],
            color="gray",
            alpha=0.2,
            label=f"Bin {i + 1}" if i == 0 else "",
        )

    ax_mu.set_xlabel(r"$\mu = \cos\theta$", fontsize=12)
    ax_mu.set_ylabel("Normalized Intensity", fontsize=12)
    ax_mu.set_title(title + " (vs. $\mu$)", fontsize=14)
    ax_mu.legend()

    # Bottom row: Intensity vs theta
    ax_theta = axes[1, idx]
    ax_theta.scatter(theta, intensity, marker="o", color="blue", label="Intensity")

    ax_theta.set_xlabel(r"$\theta$ (degrees)", fontsize=12)
    ax_theta.set_ylabel("Normalized Intensity", fontsize=12)
    ax_theta.set_title(title + " (vs. $\\theta$)", fontsize=14)
    ax_theta.legend()

# Adjust layout and save the plot
plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make room for main title
plt.savefig("intensity_vs_mu_and_theta.png", dpi=300, bbox_inches="tight")

# Show the plot
plt.show()

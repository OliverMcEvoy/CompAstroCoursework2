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

# --- First Plot: Photon Scattering Only ---
fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.scatter(
    photon_data["mu"],
    photon_data["intensity"],
    color="purple",
    marker="o",
)

# Shade the bin regions
for i in range(len(photon_data)):
    ax1.axvspan(
        photon_data["bin_edge_lower"][i],
        photon_data["bin_edge_upper"][i],
        color="gray",
        alpha=0.2,
    )

ax1.set_xlabel(r"Scattering Angle $\theta$ in Radians", fontsize=12)
ax1.set_ylabel("Normalised Intensity", fontsize=12)

plt.savefig("photon_scattering.png", dpi=300, bbox_inches="tight")
plt.show()

# --- Second Plot: Comparison of Rayleigh Blue and Other ---
fig, ax2 = plt.subplots(figsize=(8, 6))

# Shade the bin regions
for i in range(len(photon_data)):
    ax2.axvspan(
        photon_data["bin_edge_lower"][i],
        photon_data["bin_edge_upper"][i],
        color="gray",
        alpha=0.2,
    )

ax2.scatter(
    rayleigh_blue_data["mu"],
    rayleigh_blue_data["intensity"],
    color="blue",
    marker="o",
    label="Blue light",
)
ax2.scatter(
    rayleigh_other_data["mu"],
    rayleigh_other_data["intensity"],
    color="green",
    marker="^",
    label="Other colours",
)

ax2.set_xlabel(r"Scattering Angle $\theta$ in Radians", fontsize=12)
ax2.set_ylabel("Normalised Intensity", fontsize=12)
ax2.legend()

plt.savefig("rayleigh_comparison.png", dpi=300, bbox_inches="tight")
plt.show()

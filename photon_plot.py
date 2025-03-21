import pandas as pd
import matplotlib.pyplot as plt

# Load data into separate DataFrames
photons_df = pd.read_csv("photon.csv", names=["y", "theta", "iterations"])
rayleigh_other_df = pd.read_csv(
    "rayleigh_other.csv", names=["y", "theta", "iterations"]
)
rayleigh_blue_df = pd.read_csv("rayleigh_blue.csv", names=["y", "theta", "iterations"])

# Create a list of DataFrames and their corresponding titles
dataframes = [photons_df, rayleigh_blue_df, rayleigh_other_df]
titles = ["Photons", "Blue light", "Other Colours"]

# Set up the figure with subplots for each dataset
fig, axes = plt.subplots(
    3, 4, figsize=(20, 15)
)  # 3 rows (one for each dataset), 4 columns (one for each plot)
fig.suptitle("Photon and Rayleigh Scattering Analysis", fontsize=16)

# Iterate over each DataFrame and create its subplots
for idx, (df, title) in enumerate(zip(dataframes, titles)):
    # Plot 1: Density function of theta
    axes[idx, 0].hist(df["theta"], bins=30, color="blue", alpha=0.6)
    axes[idx, 0].set_xlabel("Theta (radians)")
    axes[idx, 0].set_ylabel("Density")

    # Plot 2: Histogram of iterations per photon
    axes[idx, 1].hist(df["iterations"], bins=30, color="green", alpha=0.6)
    # axes[idx, 1].set_title(f"{title}: Distribution of Iterations per Photon")
    axes[idx, 1].set_xlabel("Iterations")
    axes[idx, 1].set_ylabel("Count")

    # Plot 3: Histogram of y values per photon
    axes[idx, 2].scatter(df["y"], df["theta"], color="orange", alpha=0.7)
    # axes[idx, 2].set_title(f"{title}: Distribution of y values per Photon")
    axes[idx, 2].set_xlabel("y")
    axes[idx, 2].set_ylabel("Count")
    axes[idx, 2].set_xlim(-400, 400)

    # Plot 4: Scatter plot of y vs theta
    axes[idx, 3].scatter(df["y"], df["theta"], color="purple", alpha=0.7)
    # axes[idx, 3].set_title(f"{title}: Intensity vs y")
    axes[idx, 3].set_xlabel("y")
    axes[idx, 3].set_ylabel("Theta (radians)")
    axes[idx, 3].set_xlim(-1500, 1500)


# Adjust layout
plt.show()

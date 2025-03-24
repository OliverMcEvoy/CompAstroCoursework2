import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


def plot_binned_results_comparison(sampled_files, expected_file=None):
    """Plot multiple binned results in subplots with optional expected density overlay."""

    # Create figure with subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Load expected data once if provided
    expected_data = None
    if expected_file and os.path.exists(expected_file):
        try:
            expected_df = pd.read_csv(
                expected_file,
                header=None,
                names=["bin_start", "bin_end", "expected_density"],
            )
            expected_centers = (expected_df["bin_start"] + expected_df["bin_end"]) / 2
            expected_data = (expected_centers, expected_df["expected_density"])
        except Exception as e:
            print(f"Warning: Could not load expected density file: {e}")

    # Plot each sampled file in its own subplot
    for i, (ax, sampled_file) in enumerate(zip(axes, sampled_files)):
        # Check and load sampled data
        if not os.path.exists(sampled_file):
            raise FileNotFoundError(f"Could not find results file: {sampled_file}")

        try:
            sampled_df = pd.read_csv(
                sampled_file,
                header=None,
                names=["bin_start", "bin_end", "normalized_value"],
            )
        except Exception as e:
            raise ValueError(f"Error reading sampled file {sampled_file}: {e}")

        # Calculate bin centers and widths
        bin_centers = (sampled_df["bin_start"] + sampled_df["bin_end"]) / 2
        bin_widths = sampled_df["bin_end"] - sampled_df["bin_start"]

        # Plot sampled data as bars
        bars = ax.bar(
            bin_centers,
            sampled_df["normalized_value"],
            width=bin_widths * 0.9,  # 90% width for spacing
            align="center",
            edgecolor="navy",
            linewidth=0.5,
            color="skyblue",
            alpha=0.7,
            label="Sampled distribution",
        )

        # Plot expected density if available
        if expected_data:
            ax.plot(
                expected_data[0],
                expected_data[1],
                color="black",
                linestyle="--",
                linewidth=2,
                label="Expected density",
            )

        # Customize subplot
        ax.set_xlabel(r"$\mu$ value", fontsize=12)
        if i == 0:  # Only show y-label for leftmost plot
            ax.set_ylabel("Probability Density", fontsize=12)

        # Add sample size to title
        sample_size = sampled_file.split("_")[-1].split(".")[
            0
        ]  # Extract number from filename
        ax.set_title(f"{sample_size} samples", fontsize=12)

        # Adjust x-axis limits
        padding = bin_widths.iloc[0] * 0.5
        ax.set_xlim(
            sampled_df["bin_start"].min() - padding,
            sampled_df["bin_end"].max() + padding,
        )

        # Only show legend for first subplot
        if i == 2:
            ax.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig("Q1b.png")
    plt.show()


# Usage
plot_binned_results_comparison(
    sampled_files=[
        "direct_mapping_binned_results_2500.csv",
        "direct_mapping_binned_results_10000.csv",
        "direct_mapping_binned_results_50000.csv",
    ],
    expected_file="expected_density_bins.csv",
)

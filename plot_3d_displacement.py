import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Load CSV file
filename = "displacement_direction_lookup.csv"

df = pd.read_csv(filename, header=None, names=["mu", "z", "y", "x"])

# Extract columns
mu = df["mu"]
x = df["x"]
y = df["y"]
z = df["z"]

# Create figure with subplots
fig = plt.figure(figsize=(15, 10))

# Create 3D scatter plot (top subplot)
ax1 = fig.add_subplot(211, projection="3d")
sc = ax1.scatter(x, y, z, c=mu, marker="o", s=3, cmap="coolwarm")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.set_zlabel("Z")
ax1.set_aspect("equal")

# Create 2D scatter plots (bottom row)
ax2 = fig.add_subplot(236)
ax2.scatter(x, y, c=mu, s=3, cmap="coolwarm")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_aspect("equal")

ax3 = fig.add_subplot(235)
ax3.scatter(y, z, c=mu, s=3, cmap="coolwarm")
ax3.set_xlabel("Y")
ax3.set_ylabel("Z")
ax3.set_aspect("equal")

ax4 = fig.add_subplot(234)
ax4.scatter(x, z, c=mu, s=3, cmap="coolwarm")
ax4.set_xlabel("X")
ax4.set_ylabel("Z")
ax4.set_aspect("equal")

# Adjust layout
plt.tight_layout()

# Show plot
plt.savefig("sphere_with_subplots.png")
plt.show()

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


# Create 3D scatter plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")

# Scatter plot with color based on mu
sc = ax.scatter(x, y, z, c=mu, cmap="viridis", marker="o")

# Add color bar
cbar = plt.colorbar(sc, ax=ax, shrink=0.6, aspect=10)
cbar.set_label("mu value")

# Labels and title
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("3D Scatter Plot on Unit Sphere")
ax.set_aspect("equal")

# Show plot
plt.show()

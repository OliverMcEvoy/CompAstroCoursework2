import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Load CSV file
filename = "displacement_direction_lookup.csv"
df = pd.read_csv(filename, header=None, names=["mu", "z", "y", "x"])

# Extract columns
mu = df["mu"].values
x = df["x"].values
y = df["y"].values
z = df["z"].values

# Check unit sphere condition
radii = np.sqrt(x**2 + y**2 + z**2)
print(max(radii))
print(min(radii))
sphere_check = np.allclose(radii, 1, atol=1e-5)  # Should be close to 1
print(f"All points lie on the unit sphere: {sphere_check}")

# Check if mu maps correctly to z
mu_check = np.allclose(mu, (1 + z) / 2, atol=1e-3)
print(f"mu correctly maps to z: {mu_check}")

# Create 3D scatter plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")

# Scatter plot
sc = ax.scatter(x, y, z, marker="o", s=3, c=mu, cmap="viridis")

# Labels and title
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_aspect("auto")

# Show plot
plt.colorbar(sc, label="mu value")
plt.savefig("sphere.png")
plt.show()

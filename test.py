import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm  # Import tqdm for progress bar

# Constants
zmin = 0
zmax = 200
tau = 10  # Total optical depth
albedo = 1  # Perfect scattering
num_photons = 100
num_bins = 10
bin_edges = np.linspace(0, 1, num_bins + 1)  # Bins in mu = cos(theta)
bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
intensity = np.zeros(num_bins)

# Absorption coefficient (alpha_nu)
alpha_nu = tau / (zmax - zmin)


# Function to sample a random direction isotropically
def isotropic_direction():
    phi = 2 * np.pi * np.random.rand()  # Azimuthal angle
    mu = 2 * np.random.rand() - 1  # Cosine of polar angle (isotropic)
    theta = np.arccos(mu)  # Polar angle
    return theta, phi, mu


# Monte Carlo simulation
for _ in tqdm(range(num_photons), desc="Simulating photons", unit="photon"):
    # Initialize photon at origin
    x, y, z = 0, 0, zmin
    theta, phi, mu = isotropic_direction()  # Initial direction

    while True:
        # Calculate distance to next scattering event
        tau_scatter = -np.log(np.random.rand())  # Exponential distribution
        dz = tau_scatter / alpha_nu  # Vertical distance

        # Update position
        dx = dz * np.sin(theta) * np.cos(phi)
        dy = dz * np.sin(theta) * np.sin(phi)
        dz = dz * np.cos(theta)
        x += dx
        y += dy
        z += dz

        # Check if photon escapes
        if z >= zmax:
            # Bin the escape direction
            mu_escape = np.cos(theta)
            bin_index = np.searchsorted(bin_edges, mu_escape) - 1
            if bin_index >= 0 and bin_index < num_bins:
                intensity[bin_index] += 1
            break

        # Scattering event (isotropic)
        theta, phi, mu = isotropic_direction()

# Normalize intensity
intensity /= num_photons

# Print intensity bins
print("Intensity bins (mu, intensity):")
for mu, I in zip(bin_midpoints, intensity):
    print(f"{mu:.2f}: {I:.6f}")

# Plot intensity vs mu
plt.figure(figsize=(8, 6))
plt.plot(bin_midpoints, intensity, "o-", label="Intensity")
plt.xlabel(r"$\mu = \cos(\theta)$")
plt.ylabel("Normalized Intensity")
plt.title("Intensity as a Function of Angle")
plt.grid(True)
plt.legend()
plt.show()

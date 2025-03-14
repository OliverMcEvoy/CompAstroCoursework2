import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Get data.
data = pd.read_csv("efficient_results.csv", names=["mu", "y"])
mu_values = data["mu"]
y_values = data["y"]

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Plot all results in the first subplot.
axs[0].hist(
    mu_values, bins=300, alpha=1, label="All Results", color="orange", density=True
)
axs[0].set_title("All Results")
axs[0].set_xlabel("Mu Values")
axs[0].set_ylabel("Probability Density")
axs[0].legend()

# Plot the actual distribution.
mu_range = np.linspace(-1, 1, 1000)
theoretical_distribution = (3 / 8) * (1 + mu_range**2)
theoretical_distribution /= np.trapz(theoretical_distribution, mu_range)

axs[0].plot(mu_range, theoretical_distribution, "k--", label="Ideal Distribution")
axs[0].legend(loc="upper right")

# Scatter plot of mu vs y in the second subplot.
axs[1].scatter(
    mu_values, y_values, alpha=1, color="black", s=1, label="results of y against mu"
)

# Plot equation of mu against y for comparison.
mu_range = np.linspace(min(mu_values), max(mu_values), 100)
y_equation = (mu_range**3 + 3 * mu_range + 4) / 8
axs[1].scatter(
    mu_range, y_equation, color="red", s=10, alpha=0.5, marker="^", label="Equation Fit"
)

axs[1].set_title("Mu vs Y")
axs[1].set_xlabel("Mu Values")
axs[1].set_ylabel("Y Values")
axs[1].legend()

plt.tight_layout()
plt.show()

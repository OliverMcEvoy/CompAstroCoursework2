import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Get data.
data = pd.read_csv("results.csv")
values = data.iloc[:, 0]
labels = ["Rejected", "Accepted"]
accepted_values = values[data["accepted"] == "accepted"].values
not_accepted_values = values[data["accepted"] == "rejected"].values

fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# Plot all results in the first subplot.
axs[0].hist(
    values, bins=300, alpha=1, label="All Results", color="orange", density=True
)
axs[0].set_title("All Results")
axs[0].set_xlabel("Values")
axs[0].set_ylabel("Probability Density")
axs[0].legend()

# Plot accepted results in the second subplot.
axs[1].hist(
    accepted_values, bins=300, alpha=1, label="Accepted", color="blue", density=True
)
axs[1].set_title("Accepted Results")
axs[1].set_xlabel("Values")
axs[1].set_ylabel("Probability Density")

# plot the actual distriubiton.
mu_range = np.linspace(-1, 1, 1000)
theoretical_distribution = (3 / 8) * (1 + mu_range**2)
theoretical_distribution /= np.trapz(theoretical_distribution, mu_range)

axs[0].plot(
    mu_range,
    theoretical_distribution,
    "k--",
    label="Ideal distribution",
)
axs[1].plot(
    mu_range,
    theoretical_distribution,
    "k--",
    label="ideal distribution",
)

# Add legends to both subplots
axs[0].legend(loc="upper right")
axs[1].legend(loc="upper right")

plt.tight_layout()
plt.show()

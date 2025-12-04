import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
This program uses a monte carlo simulation to explore how variability in thermodynamics and binding free energy 
affects protein–ligand binding. By sampling experimentally yielded Kds and converting them into binding free energies,
this simulation provides insights on expected binding outcomes. This approach is especially useful for 
under-researched protein–ligand complexes, where experimental data is sparse and binding behavior is uncertain.
"""

# Load in kD data
df = pd.read_csv("gefitnib_egfr.csv")
Kd = df["Standard Value"].values
Kd *= 10**(-9)  # convert from nM -> M

# Convert kDs -> delta Gs
R = 8.31446     # gas constant (J / (mol*K))
T = 298         # temperature (Kelvin, 25 C)
deltaG = R * T * np.log(Kd)

def monte_carlo_sim(num_trials = 10000, ligand_concentration_nM = 10):
    
     # Convert ligand concentration from nM -> M
    L = ligand_concentration_nM * 10**(-9)
    
    # Randomly sample delta G values - normal distribution
    mu = deltaG.mean()
    sigma = deltaG.std()
    sampled_deltaG = np.random.normal(mu, sigma, size=num_trials)

    # Convert sampled delta G back into Kd
    sampled_Kd = np.exp(sampled_deltaG / (R * T))

    # Compute binding probability
    p_binding = L / (L + sampled_Kd)

    # Simulate binding events
    binding_events = np.random.binomial(1, p_binding)

    return sampled_deltaG, sampled_Kd, p_binding, binding_events

N = 10000   # number of trials
L_nM = 10   # ligand concentration (nM)

# Run simulation
dG, Kd, p, binding_events = monte_carlo_sim(num_trials=N, ligand_concentration_nM=L_nM)

# Plot delta G and Kd distributions
plt.hist(dG, bins=30, edgecolor='blue')
plt.title("Histogram of Sampled Binding Free Energies (ΔG)")
plt.xlabel("ΔG (kcal/mol)")
plt.ylabel("Count")
plt.savefig("FreeEnergyDistribution.png", dpi=300, bbox_inches='tight')
plt.show()

plt.hist(p, bins=30, edgecolor='blue')
plt.title("Histogram of Binding Probabilities")
plt.xlabel("Binding Probability")
plt.ylabel("Count")
plt.savefig("BindingProbabilityDistribution.png", dpi=300, bbox_inches='tight')
plt.show()

# Print simulation results
print("Monte Carlo Simulation of EGFR–Gefitinib Binding")
print(f"Trials: {N}")
print(f"Ligand concentration: {L_nM} nM")
print(f"Mean delta G: {np.mean(dG):.2f} J/mol")
print(f"Mean Kd: {np.mean(Kd)*(10**(-9)):.2f} nM")
print(f"Mean binding probability: {np.mean(p)*100:.1f}%")
print(f"Binding events: {np.sum(binding_events)} / {N}")


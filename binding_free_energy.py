import numpy as np
import pandas as pd

# load in kD data
df = pd.read_csv("gefitnib_egfr.csv")
Kd = df["Standard Value"].values
Kd *= 1e-9

# convert kDs -> delta Gs
R = 1.987e-3
T = 298
deltaG = R * T * np.log(Kd)

def monte_carlo_sim(num_trials = 10000, ligand_concentration_nM = 10):
    
     # Convert ligand concentration from nM -> M
    L = ligand_concentration_nM * 1e-9
    
    # Randomly sample delta G values - normal distribution
    mu = deltaG.mean()
    sigma = deltaG.std()
    sampled_deltaG = np.random.normal(mu, sigma, size=num_trials)
    # sampled_deltaG = np.random.choice(deltaG, size=num_trials, replace=True)

    # Convert sampled delta G back into Kd
    sampled_Kd = np.exp(sampled_deltaG / (R * T))

    # Compute binding probability
    p_binding = L / (L + sampled_Kd)

    # Generate Bernoulli binding events
    binding_events = np.random.binomial(1, p_binding)

    return sampled_deltaG, sampled_Kd, p_binding, binding_events

N = 10000   # number of trials
L_nM = 10   # ligand concentration (nM)

dG, Kd, p, binding_events = monte_carlo_sim(num_trials=N, ligand_concentration_nM=L_nM)

print("Monte Carlo Simulation of EGFR–Gefitinib Binding")
print(f"Trials: {N}")
print(f"Ligand concentration: {L_nM} nM")
print(f"Mean delta G: {np.mean(dG):.2f} kcal/mol")
print(f"Mean Kd: {np.mean(Kd)*1e9:.2f} nM")
print(f"Mean binding probability: {np.mean(p)*100:.1f}%")
print(f"Binding events: {np.sum(binding_events)} / {N}")
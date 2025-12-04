**Modeling Protein-Ligand Binding with Monte Carlo Simulations**  
**By: Jonas Kazimli, Sibhi Sakthivel, Amala Vellappillil Biju, Leonard Wei**

# Program Overview

This repository contains the code and data used for the CHEM 274B micropresentation on using Monte Carlo simulations to simulate the binding outcomes of the EGFR-Gefitinib protein-ligand complex. The simulation script `binding_free_energy.py` uses experimentally yielded Kd values to randomly sample binding free energy values to determine probability of binding and count how many binding events occur out of a specified number of trials.

# Setup Instructions

1) Clone or download this repository.
2) Install dependencies:
```bash
pip install -r requirements.txt
```

3) Run the simulation script:
```bash
python binding_free_energy.py
```

# Expected Output

Running `binding_free_energy.py` generates:
- A histogram of sampled binding free energies (`FreeEnergyDistribution.png`)
- A histogram of binding probabilities (`BindingProbabilityDistribution.png`)
- Console output summarizing:
    - number of Monte Carlo trials
    - number of binding events
    - estimated binding probability

# Dataset

The file `gefitinib_egfr.csv` contains experimentally measured dissociation constants (Kd) for the EGFR–Gefitinib complex obtained from chEMBL, a public database containing a variety of data related to how small molecules interact with their protein targets. These values serve as the foundation for sampling delta G values in the Monte Carlo simulation.

https://www.ebi.ac.uk/chembl/explore/activities/STATE_ID:S_ELsqXybdbPWYqKDsPEIA%3D%3D

# Modifying Simulation Parameters

You can adjust the number of Monte Carlo trials or modify simulation parameters by editing 
the constants in `binding_free_energy.py`:

```python
N = 10000           # number of Monte Carlo iterations
R = 8.31446         # gas constant (J/(mol·K))
T = 298             # temperature (K)
L_nM = 10           # ligand concentration (nM)
```

# Repository Contents
binding_free_energy.py              # Main simulation script
gefitinib_egfr.csv                  # Sample dataset with experimentally yielded Kds
FreeEnergyDistribution.png          # Free energy distribution plot
BindingProbabilityDistribution.png  # Binding probability distribution plot
requirements.txt                    # Python dependencies
README.md                           # Overview & setup instructions

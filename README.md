**Modeling Protein-Ligand Binding with Monte Carlo Simulations**  
**By: Jonas Kazimli, Sibhi Sakthivel, Amala Vellappillil Biju, Leonard Wei**

# Program Overview

This repository contains the code and data used for the CHEM 274B micropresentation on using monte carlo simulations to simulate the binding outcomes of the EGFR-Gefitinib protein-ligand complex. The simulation script `binding_free_energy.py` uses experimentally yielded Kd values to randomly sample binding free energy values to determine probability of binding and count how many binding events occur out of a specified number of trials.

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

# Repository Contents
binding_free_energy.py              # Main simulation script
gefitinib_egfr.csv                  # Sample dataset with experimentally yielded Kds
FreeEnergyDistribution.png          # Free energy distribution plot
BindingProbabilityDistribution.png  # Binding probability distribution plot
requirements.txt                    # Python dependencies
README.md                           # Overview & setup instructions

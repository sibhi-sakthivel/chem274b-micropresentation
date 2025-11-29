import numpy as np

# we need to choose a protein and ligand to run the simulation on

# define kD - we'll get this value from BindingDB
kd = 0

# we need to get a dataset with raw ligand concentrations
# then choose a uniform or normal distribution to randomize [L] in the simulation

# define [L] upper and lower bounds - uniform distribution
# upper 0
# lower = 0

# define [L] mean, variance - normal distribution
mean = 0
variance = 0

# define number of iterations
num_trials = 1000

# track binding events during simulation
bound = 0

# normal distribution - [L]
concentrations = np.random.normal(loc=mean, scale=variance, size=num_trials)

# or uniform distribution
# concentrations = np.random.uniform(upper, lower, num_trials)

# calculate probability of binding for each random concentration
binding_probabilities = concentrations / (concentrations + kd)

# count number of binding events
for i in range(num_trials):
    if np.random.uniform(0, 1) < binding_probabilities[i]:
        bound += 1
        
fraction_bound = bound / num_trials

print(f"{bound} out of {num_trials} proteins formed a protein-ligand complex: {fraction_bound*100}%")
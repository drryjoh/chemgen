import numpy as np 
import matplotlib.pyplot as plt

import numpy as np
import cantera as ct

def safe_divide(a, b):
    """
    Safely divides two inputs, handling cases where division by zero or near-zero values might occur.
    Works for both scalars and NumPy arrays.
    
    Args:
    - a: numerator (scalar or array)
    - b: denominator (scalar or array)
    
    Returns:
    - Division result with safe handling for near-zero denominators.
    """
    if isinstance(a, (np.ndarray, list)) or isinstance(b, (np.ndarray, list)):
        # Convert inputs to NumPy arrays for consistent handling
        a = np.array(a, dtype=np.float64)
        b = np.array(b, dtype=np.float64)
        
        # Avoid division by very small or zero values
        safe_b = np.where(np.abs(b) < 1e-20, np.nan, b)
        return np.divide(a, safe_b, out=np.zeros_like(a), where=~np.isnan(safe_b))
    else:
        # Handle scalar inputs
        if np.abs(b) < 1e-20:
            return 0.0
        else:
            return a / b

# Load the CSV data
state_files = ["bad_state_4745.npy"]

mech_file = "FFCM2_model.yaml"
gases = []
for point in range(len(state_files)):
    bad_state = np.load(state_files[point])
    gas = ct.Solution(mech_file)
    temperature_b = bad_state[0]
    pressure_b = ct.gas_constant * temperature_b * np.sum(bad_state[1:])
    X_b= bad_state[1:]/np.sum(bad_state[1:])
    gas.TPX =  temperature_b, pressure_b, X_b
    gases.append(gas)

data = np.loadtxt("l2_norm_results.csv", skiprows=1, delimiter=',')

# Extract columns
temperatures = data[:, 0]  # First column: temperatures
errors = data[:, 1:-1]       # Remaining columns: chemgen and cantera errors
print(errors.shape)
# Separate chemgen and cantera errors
chemgen = errors[:, ::2]  # Every other column starting at 0
cantera = errors[:, 1::2]  # Every other column starting at 1
print(chemgen)
# Validate shapes
#pt, reaction
print("Temperatures shape:", temperatures.shape)
print("Chemgen shape:", chemgen.shape)
print("Cantera shape:", cantera.shape)

# Plot the scatter plot
reactions_to_explore = []
errors_to_explore = []
color_idx = np.linspace(0, 1, chemgen.shape[0])
random_shift = [np.random.uniform(100) for i in range(chemgen.shape[0])]
for i in range(chemgen.shape[0]): #pt
    for k in range(chemgen.shape[1]): #reactionc
        # Calculate relative error with safe division
        error = 0
        if cantera[i, k] != 0:
            error = np.abs((chemgen[i, k] - cantera[i, k]) / cantera[i, k])
            #print(error)
            if error>0.1:
                plt.plot(temperatures[i]+random_shift[i], error, 'o', markersize=1, color=plt.cm.jet(color_idx[i]))
                reactions_to_explore.append(k)
                errors_to_explore.append(error)
                plt.text(temperatures[i]+random_shift[i], error, f"{k+1}", fontsize=8, ha='left')
                C = gases[i].X
                print(f"{gases[i].reaction(k)}")
                print(f"{gases[i].reaction(k).reaction_type}")
                #print(f" chemgen value: {chemgen[i, k]}, cantera value: {cantera[i, k]} for reaction {k+1} recalculation: {gases[i].net_rates_of_progress[k]}")
                #for reactant_species in gases[i].reaction(k).reactants:
                #    idx = gases[i].species_names.index(reactant_species)
                #    print(f"reactant: {reactant_species}: {C[idx]}")
                #for product_species in gases[i].reaction(k).products:
                #    idx = gases[i].species_names.index(product_species)
                #    print(f"product: {product_species}: {C[idx]}")
        else:
            error = 0  # Handle zero denominator safely
# Combine, sort, and unzip
sorted_pairs = sorted(zip(errors_to_explore, reactions_to_explore), key=lambda x: x[0], reverse=True)
sorted_errors, sorted_things = zip(*sorted_pairs)

# Convert to lists (optional, as zip returns tuples)
sorted_errors = list(sorted_errors)
sorted_things = list(sorted_things)

# Output the results
#print("Sorted Errors:", sorted_errors)
print("Sorted Things:", np.unique(sorted_things))

plt.show()
## Plot the scatter points for visualization
#plt.xlabel("Reaction Index")
#plt.ylabel("Temperature")
#plt.title("Annotated Errors")
#plt.show()


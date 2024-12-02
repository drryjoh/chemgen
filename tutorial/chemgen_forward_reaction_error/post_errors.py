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
state_files = ["bad_state_1662.npy", "bad_state_4577.npy", "bad_state_4831.npy", "bad_state_7395.npy"]
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
    print(pressure_b/101325.0)

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
for i in range(chemgen.shape[0]): #pt
    for k in range(chemgen.shape[1]): #reaction
        # Calculate relative error with safe division
        error = 0
        if cantera[i, k] != 0:
            error = np.abs((chemgen[i, k] - cantera[i, k]) / cantera[i, k])
            #print(error)
            plt.plot(temperatures[i], error, 'ok', markersize=1)
            if error>0.2:
                plt.text(temperatures[i], error, f"{k+1}", fontsize=8, ha='left')
                #print(gases[i].forward_rate_constants[k])
                C = gases[i].X
                print(f" chemgen value: {chemgen[i, k]}, cantera value: {cantera[i, k]} for reaction {k+1}")
                print(f"reaction information:\n {gases[i].reaction(k)}")
                #for reactant in gases[i].reaction(k).reactants:
                #    print(C[gases[i].species_names.index(reactant)])
                #for product in gases[i].reaction(k).products:
                #    print(C[gases[i].species_names.index(product)])
        else:
            error = 0  # Handle zero denominator safely
plt.show()
## Plot the scatter points for visualization
#plt.xlabel("Reaction Index")
#plt.ylabel("Temperature")
#plt.title("Annotated Errors")
#plt.show()


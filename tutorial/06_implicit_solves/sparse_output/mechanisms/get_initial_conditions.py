import cantera as ct

# Temperature: 700-1600 K
# Pressure: 1-2 atm

flag = 2

if flag == 0:
    mech_file = "../ffcm2_h2.yaml"
    fuel_name = "H2"
elif flag == 1:
    mech_file = "../FFCM2_model.yaml"
    fuel_name = "C2H4"
    # Jacobian: # nonzero entries = 3256, # zero entries = 6153, total # entries = 9409
    # Jacobian: sparsity percentage = 65.3948%
    # LU decomposition of permuted Jacobian: # nonzero entries = 3965, # zero entries = 5444, total # entries = 9409
    # LU decomposition of permuted Jacobian: sparsity percentage = 57.8595%
elif flag == 2:
    mech_file = "3-methylheptane-c8h18-3-1378-8143.yaml"
    fuel_name = "c8h18-3"
    # Jacobian: # nonzero entries = 26216, # zero entries = 1875425, total # entries = 1901641
    # Jacobian: sparsity percentage = 98.6214%
    # LU decomposition of permuted Jacobian: # nonzero entries = 50017, # zero entries = 1851624, total # entries = 1901641
    # LU decomposition of permuted Jacobian: sparsity percentage = 97.3698%
elif flag == 3:
    mech_file = "hydrogen-10-28.yaml"
    fuel_name = "H2"
elif flag == 4:
    mech_file = "jetA-detailed-NOx-203-1589.yaml"
    fuel_name = "POSF10325"
elif flag == 5:
    mech_file = "butane-230-2461.yaml"
    fuel_name = "C4H10"
elif flag == 6:
    mech_file = "ic8-874-6864.yaml"
    fuel_name = "IC8H18"
elif flag == 7:
    mech_file = "md5d-2649-10487.yaml"
    fuel_name = "MD5D"
elif flag == 8:
    mech_file = "md-nc7-3787-10264.yaml"
    fuel_name = "nc7h16"
elif flag == 9:
    mech_file = "mmc5-7171-38324.yaml"
    fuel_name = "c20h42-2"


else:
    raise ValueError


gas = ct.Solution(mech_file)

air = "O2:0.21,N2:0.79"
# air = "O2:1"

fuel = fuel_name + ":1"

gas.set_equivalence_ratio(phi=1.0, fuel=fuel, oxidizer=air)


def print_mole_fraction(species_name):
    print("""    - name: {0}
      value: {1}""".format(species_name, gas[species_name].X[0]))

# Print
print("  MoleFraction:")
print_mole_fraction("O2")
print_mole_fraction("N2")
print_mole_fraction(fuel_name)

# breakpoint()

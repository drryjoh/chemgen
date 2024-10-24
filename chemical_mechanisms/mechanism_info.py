import cantera as ct

# Load the mechanism file
# Replace 'gri30.xml' with the path to your mechanism file if necessary
gas = ct.Solution('FFCM2_model.yaml')

# Initialize a dictionary to count different reaction types
reaction_types_count = {}

# Loop over all reactions in the mechanism
for i in range(gas.n_reactions):
    reaction = gas.reaction(i)
    reaction_type = reaction.reaction_type 

    # Increment the count of this reaction type
    if reaction_type not in reaction_types_count:
        reaction_types_count[reaction_type] = 1
    else:
        reaction_types_count[reaction_type] += 1

# Print out the reaction types and their counts
print("Reaction types and their counts:")
for reaction_type, count in reaction_types_count.items():
    print(f"{reaction_type}: {count}")

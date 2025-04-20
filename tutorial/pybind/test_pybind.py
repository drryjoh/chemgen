#!python3
import chemwrapper

species = [1.0, 2.0, 3.0, 4.0, 5.0]
scaled = chemwrapper.scale_species(species, 2.0)
print(scaled)  # Output: [2.0, 4.0, 6.0, 8.0, 10.0]


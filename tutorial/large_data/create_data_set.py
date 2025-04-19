#!python3
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

gas = ct.Solution("H2_FFCM1_O3_Ar.yaml")
data = []
for name in [ "T", "p", *gas.species_names]:
    print(f"processing {name}")
    f = open(f'0.0002/{name}','r')
    reading = False
    read_number_of_points = False
    number_of_points = 0
    points = 0
    d = []
    for line in f:
        if reading and points < number_of_points and "(" not in line:
            d.append(float(line.strip(" ").strip("\n")))
            points+=1
        if "internalField   nonuniform List<scalar>" in line:
            print("start")
            read_number_of_points = True
        elif read_number_of_points:
            number_of_points = int(line.strip(" ").strip("\n"))
            print(number_of_points)
            read_number_of_points = False
            reading = True
    d = np.array(d)
    if len(d) != number_of_points:
        print("error!")
    data.append(d)
data = np.array(data).transpose()

f = open("consolidated.csv","w")
for line in data:
    str_line = [str(d) for d in line]
    f.write("{0}\n".format(", ".join(str_line)))



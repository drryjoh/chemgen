#!python3
from cantera import ck2yaml
ck2yaml.convert(input=['ffcnoxmod.inp '], thermo='thermnox.dat', out_name='mechanism.yaml')


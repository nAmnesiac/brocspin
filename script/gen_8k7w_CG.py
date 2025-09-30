import os 
import numpy as np
from HyresBuilder import RNAbuilder

sequence = "GGGACGGUCGGGUCCAGAUAUUCGUAUCUGUCGAGUAGAGUGUGGGCUC"
filepath = "/home/allenchen/brocspin/structures/CG_pdb/8k7w_CG_do.pdb"

print(sequence)
print(filepath)

pdb_dir = os.path.dirname(filepath)
if not os.path.exists(pdb_dir):
     os.makedirs(pdb_dir)

print(f'Generating CG PDB: {sequence} at path {filepath}')
RNAbuilder.build(sequence, filepath)
print('Done')

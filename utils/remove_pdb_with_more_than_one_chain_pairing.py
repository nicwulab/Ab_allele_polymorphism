# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

from Bio.Seq import Seq
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from Bio.PDB.DSSP import DSSP
from Bio import SeqIO
from Bio.SeqUtils import seq1

import pandas as pd
import os

# Read TSV file into DataFrame
df = pd.read_table('pdb_with_only_one_chain_pairing.tsv')

pdb_list = df['pdb'].tolist()
print(pdb_list)
print(len(pdb_list))
set_pdb_list = set(pdb_list)
deleted = list()
folder = '/Users/natalieso/Downloads/20230217_0084705/'
for root, dirs, files in os.walk(folder):
    for name in files:
        if name.endswith((".pdb")):
            pdb_name = name.split('.')[0]
            if pdb_name not in set_pdb_list:
                os.remove(os.path.join(root, name))
                deleted.append(pdb_name)

file = open('deleted_items.txt', 'w')
for d in deleted:
    file.write(d + "\n")
file.close()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

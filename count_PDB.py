#!/usr/bin/python3
import os
import sys
import pandas as pd

PDB_all = eval(open('results/pdb_to_paratope/dssp_pdb_seq_location.txt', 'r').read())
PDB_all = set(list(PDB_all.keys()))
PDB_nan_antigen = eval(open('results/nan_antigen_chain.txt', 'r').read())
PDB_nan_antigen = set(PDB_nan_antigen)
PDB_yes_antigen = PDB_all - PDB_nan_antigen
PDB_final = set(pd.read_csv('results/allele_var_info_table_final.csv')['pdb_id'])
print ('Total PDB: %i' % len(PDB_all))
print ('PDB with none antigen: %i' % len(PDB_nan_antigen))
print ('PDB with antigen: %i' % len(PDB_yes_antigen))
print ('PDB final (original): %i' % len(PDB_final))
print ('PDB final (with antigen): %i' % len(PDB_final-PDB_nan_antigen))

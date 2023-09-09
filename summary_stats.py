#/usr/bin/python3
import sys
import pandas as pd
from collections import Counter

def loc_classification(loc):
  try:
    loc = int(loc)
  except:
    loc = int(loc[0:-1])
  if loc >= 31 and loc <= 35:
    return ('CDRH1')
  elif loc >= 50 and loc <= 65:
    return ('CDRH2')
  else: 
    return ('FW')

def main():
  df = pd.read_csv('results/allele_var_info_table_final.csv')
  PDBs  = list(df['pdb_id'])
  print ("Total PDB:", len(set(PDBs)))
  df = df[df['DDG_binding']>0]
  PDBs  = list(df['pdb_id'])
  genes = list(df['gene'])
  locations = list(df['location'])
  aas   = list(df['amino_acid_original'])
  muts  = ['-'.join(mut) for mut in zip(genes, aas, locations)]
  mut_class = list(map(loc_classification, locations))
  print ('ddG > 0:')
  print ("# of PDB:", len(set(PDBs)))
  print ('# of genes:', len(genes))
  print ('# of unique polymorphism:', len(set(muts)))
  print ("# of unique genes:", len(set(genes)))
  print ("Gene family:", ",".join(sorted(set(map(lambda x:x.rsplit('-')[0], genes)))))
  print (Counter(mut_class))

if __name__ == "__main__":
  main()


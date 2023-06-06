#!/usr/bin/python3
import sys
import numpy as np
import pandas as pd

def excel_to_dict(filename, index_colname):
  df = pd.read_excel(filename, sheet_name="Sheet1",index_col=index_colname)
  return df.to_dict()
  
def antigen_ID_mapping(df, species_to_ID, antigen_to_ID):
  ID = ''
  try:
    ID = antigen_to_ID['ID'][df['antigen_name']]
  except:
    ID = species_to_ID['ID'][df['antigen_species'].rsplit('|')[0].rstrip()]
  if ID == '':
    print ('error')
    sys.exit()
  return ID

def resolution_mapping(df, PDB_to_resolution):
  resolution = ''
  try: 
    resolution = PDB_to_resolution['resolution'][df['pdb_id']]
  except:
    resolution = df['resolution']
  return resolution

def example_classification(df, known_examples):
  if df['pdb_id'] in known_examples['Mut'].keys() and \
  known_examples['Mut'][df['pdb_id']] == df['amino_acid_original']+str(df['location'])+df['variants']:
    return 'yes'
  else:
    return 'no'
  
def location_lower_case(df, colname):
  try :
    return str(int(df[colname]))
  except:
    return df[colname][0:-1]+df[colname][-1].lower()


def add_antigen_ID(df, species_to_ID, antigen_to_ID, PDB_to_resolution, abbreviate_ID, known_examples):
  df['antigen_ID'] = df.apply(antigen_ID_mapping, axis=1, args=(species_to_ID, antigen_to_ID,))
  df['resolution'] = df.apply(resolution_mapping, axis=1, args=(PDB_to_resolution,))
  df['known'] = df.apply(example_classification, axis=1, args=(known_examples,))
  df['location'] = df.apply(location_lower_case, axis=1, args=("location",))
  df['dssp_location'] = df.apply(location_lower_case, axis=1, args=("dssp_location",))
  df['antigen_abbrev'] = df['antigen_ID'].map(abbreviate_ID['abbrev'])
  df['antigen_class'] =  df['antigen_ID'].map(abbreviate_ID['class'])
  df['DDG_binding'] = df['DDG'].astype('float64') - df['DDG_antibody'].astype('float64')
  df = df[df['DDG']!="-nan\n"]
  df = df.drop(['Unnamed: 0'], axis=1)
  df.to_csv('results/allele_var_info_table_final.csv',index=False)
  
def main():
  df = pd.read_csv('results/add_dssp_area.csv')
  species_to_ID = excel_to_dict('doc/species_ID.xlsx', 'species')
  antigen_to_ID = excel_to_dict('doc/antigen_name_ID.xlsx', 'antigen_name')
  PDB_to_resolution = excel_to_dict('doc/resolution.xlsx', 'PDB')
  abbreviate_ID = excel_to_dict('doc/abbreviate_ID.xlsx', 'ID')
  known_examples = excel_to_dict('doc/known_examples.xlsx', 'PDB')
  df = add_antigen_ID(df, species_to_ID, antigen_to_ID, PDB_to_resolution, abbreviate_ID, known_examples)

if __name__ == "__main__":
  main()

# Ab_allele

### Dependencies ###
* [python](https://www.python.org/) (version 3.9)
* [PyIR](https://github.com/crowelab/PyIR)
* [BioPython](https://github.com/biopython/biopython)
* [Pandas](https://pandas.pydata.org/)
* [ANARCI](https://github.com/oxpig/ANARCI)
## Dependencies Installation ##
Install dependencies by conda:

```
conda create -n Abs \
  python=3.9 \
  biopython \
  pandas \
  igblast \
  anarci
```

## Input files 

* List of pdb files and summary file from downloaded from https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?ABtype=All&method=All&species=HOMO+SAPIENS&resolution=&rfactor=&antigen=All&ltype=All&constantregion=All&affinity=All&isin_covabdab=All&isin_therasabdab=All&chothiapos=&restype=ALA&field_0=Antigens&keyword_0=#downloads
scroll down to the bottom of the page to download an archived zip file
## Local PyIR setup

### PyIR: An IgBLAST wrapper and parser

``pip3 install crowelab_pyir``

Database set up in pyir library directory

``pyir setup``

### Manually install IMGT REF database

1. Sequence download from  http://www.imgt.org/vquest/refseqh.html#VQUEST

2. Copy and paste, save as fasta (save all V gene in one file; all D gene in one file; all J gene in one file)

3. Clean data (raw edit_imgt_file.pl can be found on igblast-1.17.1xxx/bin)

``edit_imgt_file.pl imgt_database/human_prot/imgt_raw/IGV.fasta > imgt_database/human_prot/IGV.fasta``

4. Create database (use "-dbtype prot" for protein sequence, use "-dbtype nucl" for DNA sequence). For example:

``makeblastdb -parse_seqids -dbtype prot -in imgt_database/human_prot/IGV.fasta``   

``makeblastdb -parse_seqids -dbtype nucl -in imgt_database/human_nuc/IGV.fasta``

5. Run PyIR for igBlast

- see [PyIR.py](./code/_PyIR_.py)

## PDB to Paratope
1. Run DSSP to find list of amino acid locations that interact with antigens  
``python pdb_to_paratope.py {pdb_dir} {summary_file} {mkdssp_dir}``    
for example: 
``python pdb_to_paratope.py /Users/natalieso/Downloads/20230217_0084705/ data/20230217_0084705_summary.tsv /Users/natalieso/Downloads/dssp-3.1.4/mkdssp``
    - Input file:
      - List of pdb files and summary file
    - Output files:
      - [results/pdb_to_paratope/pdb_to_paratope_result_with_amino_acid_name.txt](results/pdb_to_paratope/pdb_to_paratope_result_with_amino_acid_name.txt)
      
## Allelic Variant at Paratope
1. Extract PDB sequences from PDB files  
``python utils/get_sequence_from_pdb_file.py {pdb_dir} {output_file} {summary_file}``    
for example: 
``python utils/get_sequence_from_pdb_file.py /Users/natalieso/Downloads/20230217_0084705/ results/pdb_to_paratope/pdb_sequences.fasta data/20230217_0084705_summary.tsv``
  - Input file:
      - List of pdb files and summary file
      - Output file:
        - [results/pdb_to_paratope/pdb_sequences.fasta](results/pdb_to_paratope/pdb_sequences.fasta)
2. Remove X's in fasta file  
``python utils/remove_x_in_fasta.py``
  - Input file:
    - [results/pdb_to_paratope/pdb_sequences.fasta](results/pdb_to_paratope/pdb_sequences.fasta)
  - Output file:
    - [results/pdb_to_paratope/pdb_sequences.fasta](results/pdb_to_paratope/pdb_sequences.fasta)
3. Run PyIR  
``pyir results/pdb_sequences.fasta --sequence_type prot --legacy --germlineV imgt_database/human_prot/IGV.fasta -s human``  
This generates a zip file. Please unzip the json file, rename it to {pdb_sequences.json} and move it to the results/pdb_to_paratope directory.  

4. Get list of DSSP positions for all PDB files   
``python utils/get_dssp.py {pdb_dir} {summary_file} {mkdssp_dir}``  
for example: 
``python utils/get_dssp.py /Users/natalieso/Downloads/20230217_0084705/ data/20230217_0084705_summary.tsv /Users/natalieso/Downloads/dssp-3.1.4/mkdssp``
   - Output file:
       - [results/pdb_to_paratope/dssp_pdb_seq_location.txt](results/pdb_to_paratope/dssp_pdb_seq_location.txt)
5. Get allelic variant at paratope  
``python allelic_variant_at_paratope.py`` 
   - Output file:
     - [results/allelic_variant_at_paratope/allelic_variant_at_paratope_result.csv](results/allelic_variant_at_paratope/allelic_variant_at_paratope_result.csv)
     - [results/allelic_variant_at_paratope/allelic_variant_at_paratope_result_with_dssp_loc_df.csv](results/allelic_variant_at_paratope/allelic_variant_at_paratope_result_with_dssp_loc_df.csv)
     - [results/allelic_variant_at_paratope/allelic_variant_at_paratope_result_targets_df.csv](results/allelic_variant_at_paratope/allelic_variant_at_paratope_result_targets_df.csv)
     - [results/allelic_variant_at_paratope/allelic_variant_at_paratope_result_splitted_df.csv](results/allelic_variant_at_paratope/allelic_variant_at_paratope_result_splitted_df.csv)
     
## Compute DDG

1. Run FoldX on list of PDB. Run this in the directory where you want to store the FoldX outputs.   
``python compute_ddG_foldx_script.py {pdb_dir}``   
for example:
``python compute_ddG_foldx_script.py /Users/natalieso/Downloads/20230217_0084705/``
    - Output files:
      - List of FoldX outputs containing DDG values

2. Parse FoldX output and save DDG to CSV file.   
``python compute_ddG.py {foldx output directory}``   
for example:
``python compute_ddG.py /Users/natalieso/Downloads/foldx_remote/``
    - Output file:
      - [empty_ddg_rows.csv](empty_ddg_rows.csv)
      - [FoldX.csv](FoldX.csv)  
FoldX cannot handle any positions with non-integer numbering. FoldX.csv only 
includes DDG values when the mutation target location is an integer. 
empty_ddg_rows.csv is the list where DDG was not calculated. The PDB files have 
been modified to integer-only locations, saved to /data/modified_pdb_files, and then 
run through FoldX separately. The resulting DDG values were inputted to 
results/compute_ddG/compute_ddG_result.csv
    - Output file:
      - [results/compute_ddG/compute_ddG_result.csv](results/compute_ddG/compute_ddG_result.csv)

## Epitope Identification

1. Run epitope_identification.py to get list of epitope along with a group ID for each entry in results/compute_ddG/compute_ddG_result.csv.   
``python epitope_identification.py {pdb_dir} {mkdssp_dir}``   
for example:
``python epitope_identification.py /Users/natalieso/Downloads/20230217_0084705/ /Users/natalieso/Downloads/dssp-3.1.4/mkdssp``
    - Input file:
      - [results/compute_ddG/compute_ddG_result.csv](results/compute_ddG/compute_ddG_result.csv)

    - Output files:
      - [results/epitope_idenfication/epitope_identification.csv](results/epitope_idenfication/epitope_identification.csv)
      - [results/epitope_idenfication/epitope_identification_with_group_id.csv](results/epitope_idenfication/epitope_identification_with_group_id.csv)  

## Baseline Variation

1. Run baseline_variation.py to get baseline variation with [data/anarci_igv_output.csv_H.csv](data/anarci_igv_output.csv_H.csv) and [data/anarci_igv_output.csv_KL.csv](data/anarci_igv_output.csv_KL.csv) for each entry in results/epitope_idenfication/epitope_identification_with_group_id.csv.   
``python baseline_variation.py``   
    - Input file:
      - [data/all_paired_antibodies_from_GB_v6.csv](data/all_paired_antibodies_from_GB_v6.csv)
      - [data/anarci_igv_output.csv_H.csv](data/anarci_igv_output.csv_H.csv)
      - [data/anarci_igv_output.csv_KL.csv](data/anarci_igv_output.csv_KL.csv)
      - [results/epitope_idenfication/epitope_identification_with_group_id.csv](results/epitope_idenfication/epitope_identification_with_group_id.csv)

    - Output files:
      - [results/epitope_identification_with_baseline_variation.csv](results/epitope_identification_with_baseline_variation.csv)

## Compute Antibody DDG

1. Run utils/remove_chains.py to remove chains other than antibody chains for all PDBs  
``python utils/remove_chains.py {pdb_dir} {summary_file} {output_dir}``
2. Run FoldX on list of PDB with only antibody chains. Run this in the directory where you want to store the FoldX outputs.   
``python compute_ddG_foldx_script.py {pdb_dir}``   
for example:
``python compute_ddG_foldx_script.py /Users/natalieso/Downloads/pdb_antibodies_only/``
    - Output files:
      - List of FoldX outputs containing DDG values

3. Parse FoldX output and save DDG to CSV file.   
``python compute_ddG_antidody.py {foldx output directory}``   
for example:
``python compute_ddG_antidody.py /Users/natalieso/Downloads/dir_antibodies_only_foldx_output/``
    - Output file:
      - [empty_ddg_rows.csv](empty_ddg_rows.csv)
      - [FoldX.csv](FoldX.csv)  
FoldX cannot handle any positions with non-integer numbering. FoldX.csv only 
includes DDG values when the mutation target location is an integer. 
empty_ddg_rows.csv is the list where DDG was not calculated. The PDB files have 
been modified to integer-only locations, saved to /data/modified_pdb_files_antibody_only, and then 
run through FoldX separately. The resulting DDG values were inputted to 
results/epitope_identification_with_antibody_only_ddgs.csv
    - Output file:
      - [results/epitope_identification_with_antibody_only_ddgs.csv](results/epitope_identification_with_antibody_only_ddgs.csv)

## Add DSSP RSA

1. Run add_antibody_dssp.py to get RSA from DSSP for each entry in results/epitope_identification_with_antibody_only_ddgs.csv.   
``python add_antibody_dssp.py {pdb_dir} {mkdssp_dir}``   
for example:
``python add_antibody_dssp.py /Users/natalieso/Downloads/20230217_0084705/ /Users/natalieso/Downloads/dssp-3.1.4/mkdssp``
    - Input file:
      - [results/epitope_identification_with_antibody_only_ddgs.csv](results/epitope_identification_with_antibody_only_ddgs.csv)

    - Output files:
      - [results/add_dssp_area.csv](results/add_dssp_area.csv)
      


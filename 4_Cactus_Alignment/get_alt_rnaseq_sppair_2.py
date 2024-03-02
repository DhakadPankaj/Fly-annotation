#!/usr/bin/env python

#Get altternate RNA seq species pair (species_Without_RNA_seq ALTspecies)
#Usage: python3 get_alt_rnaseq_sppair.py folder_with_cactus_files genome_info_file

import sys
import os
import pandas as pd

print(sys.argv)      
folder_with_cactus_files = sys.argv[1]  
genome_info_file = sys.argv[2]

genome_info = pd.read_csv(genome_info_file, index_col=None, sep='\t')

with open(f"{folder_with_cactus_files}/alt_rna_mapping_sp.tsv",'w') as f:
    pass

for file in os.listdir(folder_with_cactus_files):
	if file.endswith("_tree.nw"):
		cactus_file = open(f"{folder_with_cactus_files}/{file}")
		line = cactus_file.readlines()
		if len(line[5]) > 1:
        		tmp_list = line[5].replace("#Intron_bam: ","").split(",")
        		for intron_pair in tmp_list:
        			sp_without_rnaseq_hal = intron_pair.split(";")[0].rstrip("\n").strip(" ")
        			alt_sp_hal = intron_pair.split(";")[1].rstrip("\n").strip(" ")
        			sp_without_rnaseq = genome_info.loc[genome_info['hal_name'] == sp_without_rnaseq_hal, 'species_name'].item()
        			assembly = genome_info.loc[genome_info['hal_name'] == sp_without_rnaseq_hal, 'assembly'].item()
        			alt_sp = genome_info.loc[genome_info['hal_name'] == alt_sp_hal, 'species_name'].item()
        			with open(f"{folder_with_cactus_files}/alt_rna_mapping_sp.tsv", 'a') as f:
                			print(sp_without_rnaseq + "\t" + sp_without_rnaseq_hal + "\t" + assembly \
                          		+ "\t" + alt_sp+"\t" + alt_sp_hal, file=f)
  

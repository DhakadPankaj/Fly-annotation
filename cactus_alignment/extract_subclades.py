#Script for automated setup for the cactus subalignment
#it uses pandas and ete3 module for the tree readind and parsing
#Usage- python3 extract_subclades.py path_to_treefile path_to_genome_info_file

import sys
import random
import numpy as np
import pandas as pd
from ete3 import Tree, faces, AttrFace,  TreeStyle

print(sys.argv)     
treefile = sys.argv[1] 
genome_info_file = sys.argv[2]

#load tree in newick format
tree = Tree(treefile, format=1)

#defines ancestor of these species as tree outgroup
ancestor = tree.get_common_ancestor("MUSCA_DOMESTICA","GLOSSINA_MORSITANS")
tree.set_outgroup(ancestor)

#load all genome info file(file str: species_name assembly_name busco_score ref_quality cactus_name) 
genome_info = pd.read_csv(genome_info_file, index_col=None, sep='\t')

#make a list of annotated species
annotated_species = genome_info.loc[genome_info['annotated'] == "Yes", 'hal_name'].tolist()

#make a list of RNA_seq species
rnaseq_species = genome_info.loc[genome_info['RNA_seq'] == "Yes", 'hal_name'].tolist()

#List of species not included in annotation pipeline
outgroup_species = ["GLOSSINA_MORSITANS","MUSCA_DOMESTICA","PAYKULLIA_MACULATA", "POLLENIA_ANGUSTIGENA", "TACHINA_FERA",\
                   "PROTOCALLIPHORA_AZUREA", "STOMORHINA_LUNATA", "SARCOPHAGA_BULLATA", "CIRRULA_HIANS",\
                    "EPHYDRA_GRACILIS","SCATELLA_TENUICOSTA","DIASTATA_REPLETA","BRAULA_COECA"]


#add features to the tree
def tree_annotation(tree):
    count=0
    for leaf in tree.traverse():
    #add annotation feature and Ancestor name 
        if leaf.is_leaf():
            if leaf.name in annotated_species:
                leaf.add_features(annotated=True, confidence=1)
            elif leaf.name in outgroup_species:
                leaf.add_features(outgroup_sp=True, confidence=1)
            else:
                leaf.add_features(annotated=False, confidence=1)
        if leaf.is_leaf():
            if leaf.name in rnaseq_species:
                leaf.add_features(RNA_seq=True, confidence=1)
            else:
                leaf.add_features(RNA_seq=False, confidence=1)
        else:
            count = count + 1
            leaf.add_features(name=f"Anc{count}")
            
tree_annotation(tree)

print(len(tree.search_nodes(RNA_seq=True)))

#prints total number of labelled annotated species in whole tree
print ("This tree has " + str(len(tree.search_nodes(annotated=True))) + " Annotated species.")

def rev_sublist(sub_list, test_list):
    #checks if sublist is a part of list (Note: Returns False if sublist is part of list)
    if(all(x in test_list for x in sub_list)):
        return False
    else:
        return True

#set expected total number of genes in reference species
#Dmel chosen as best annotated species(newValue == total genes in Dmel)
newValue = 17868
def criterion(dict_item):
    key, value = dict_item # unpack
    return abs(value - newValue)
    
#Get the reference species for a subclade
def get_ref_sp(node, genome_info, busco_relaxation = 0.01):
    #input: node- tree/subtree and file/(dataframe/table) containing species name and busco score
    busco_dict={}
    total_genes_dict={}
    for sp in node.search_nodes(annotated=True):
        sp_name = str(sp).replace("--", "").strip("\n")
        busco_score = genome_info.loc[genome_info['hal_name'] == sp_name, 'busco_annotated'].item()           
        busco_dict[sp_name] = busco_score
        total_genes = genome_info.loc[genome_info['hal_name'] == sp_name, 'total_genes'].item()
        total_genes_dict[sp_name] = total_genes
    max_busco_sp = max(busco_dict, key=busco_dict.get)
    max_total_genes_sp = max(total_genes_dict, key=total_genes_dict.get)
    new_dict={}
    for sp, busco in busco_dict.items():
        if busco_dict[sp] >= busco_dict[max_busco_sp] - busco_relaxation:
            new_dict[sp] = total_genes_dict[sp]
    return min(new_dict.items(), key=criterion)[0]  #returns a key(ref_sp) having value closest to newValue 

#choose RNA seq of closest species for species without RNA seq
maxdist = tree.get_distance("DROSOPHILA_MELANOGASTER", target2="DROSOPHILA_SIMULANS", topology_only=False)
#print(maxdist)
def alt_rnaseq(tree, node, maxdist, ref_sp, rnaseq_species):
    nonrnaseq_sp = []
    rnaseq_sp = []
    for leaf in node.traverse():
        if leaf.is_leaf:
            if leaf.name not in rnaseq_species and leaf.name != ref_sp and (leaf.name).startswith("Anc") == False:
                nonrnaseq_sp.append(leaf.name)
            if leaf.name in rnaseq_species:
                rnaseq_sp.append(leaf.name)
    alt_rna_sp = ""
    for sp1 in nonrnaseq_sp:
        dist_dict={}
        for sp2 in rnaseq_sp:
            dist = tree.get_distance(sp1, target2=sp2, topology_only=False)
            dist_dict[sp2] = dist
        if dist_dict[min(dist_dict, key=dist_dict.get)] <= maxdist:
            alt_sp = min(dist_dict, key=dist_dict.get)
            alt_rna_sp = alt_rna_sp + sp1 + ";" + alt_sp + ", "
    return alt_rna_sp
    
            
    #leaf_names = node.get_leaf_names(is_leaf_fn=None)
    #rnaseq_sp = node.search_nodes(RNA_seq=True)
    
    
    
#Make a cactus config file containig subtree on line1, Root on line3, Ref_sp on line4, and cactus species name and location
#of their assemblies
def write_treefile_cactus(tree, node, ref_sp, genome_info):
    intron_bam_pair = alt_rnaseq(tree, node, maxdist, ref_sp, rnaseq_species)
    root_name = node.name
    other_annotated_sp = ""
    for count in range(len(node.search_nodes(annotated=True))):
        new_sps = str(node.search_nodes(annotated=True)[count]).replace("--", "").strip("\n") + ", "
        other_annotated_sp = other_annotated_sp + new_sps
    node.write(format=1, outfile=f"{root_name}_{ref_sp}_tree.nw")
    with open(f"{root_name}_{ref_sp}_tree.nw", 'a') as f:
        if len(intron_bam_pair) >= 1:
            print("\n\n#Root: " + str(node.name) +"\n#Reference_species: " + str(ref_sp) \
              +"\n#Annotated_species: " + other_annotated_sp.rstrip(", ") \
              + "\n#Intron_bam: " + intron_bam_pair.rstrip(", ") + "\n", file=f)
        else:
            print("\n\n#Root: " + str(node.name) +"\n#Reference_species: " + str(ref_sp) \
              +"\n#Annotated_species: " + other_annotated_sp.rstrip(", ") + "\n", file=f)
        for sp in node.get_leaf_names(is_leaf_fn=None):
            sp_name = str(sp).replace("--", "").strip("\n")
            assembly = genome_info.loc[genome_info['hal_name'] == sp_name, 'assembly'].item()
            cactus_name = genome_info.loc[genome_info['hal_name'] == sp_name, 'cactus_name'].item()
            print(cactus_name + " ../CAT_genomes/" + assembly, file=f)

#Walk on tree and find subclades containing ref_sp and then collapse current node.(i.e. remove all non reference species)
def collapse_tree(dist):
    prev_leaf = []
    for node in tree.iter_descendants("preorder"):
        if (len(node.search_nodes(annotated=True)) >= 1 and len(node) >= 2) \
        and rev_sublist(node.get_leaf_names(), prev_leaf)\
        and node.get_farthest_leaf(topology_only=False)[1] <= dist and len(node.search_nodes(annotated=False)) >= 1 \
        and len(node.search_nodes(outgroup_sp=True)) == 0:
            prev_leaf = node.get_leaf_names()
            if len(node.search_nodes(annotated=True)) == 1: 
                ref_sp=str(node.search_nodes(annotated=True)[0]).replace("--", "").strip("\n")
                write_treefile_cactus(tree, node, ref_sp, genome_info)                
                #remove already included non_ref species and collapse subtree(node) to ref species
                for sp in node.get_leaf_names(is_leaf_fn=None):
                    if sp != ref_sp:
                        sp = tree.search_nodes(name=sp)[0]
                        sp.delete()
            elif len(node.search_nodes(annotated=True)) > 1:
                #select species with highest busco as a reference
                root_name = node.name
                ref_sp = get_ref_sp(node, genome_info)
                write_treefile_cactus(tree, node, ref_sp, genome_info)
                for sp in node.get_leaf_names(is_leaf_fn=None):
                    if sp != ref_sp:
                        sp = tree.search_nodes(name=sp)[0]
                        sp.delete()

#Main function
def get_annot_subclade(tree, min_dist = 0.15, max_dist = 0.35, step_size = 0.05):
    #get the subclades containing atleast one annotated species
    for dist in np.arange(min_dist, max_dist, step_size):
        #dist = maximum distance between root of a subclade and farthest leaf in that subclade
        collapse_tree(dist)
    print(tree.get_ascii(show_internal=True))

#call to main function
get_annot_subclade(tree, 0.15, 0.45)

#print(t_draw.get_ascii(show_internal=True))

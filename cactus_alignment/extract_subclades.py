#Script for automated setup for the cactus subalignment
#it uses pandas and ete3 module for the tree readind and parsing

import random
import numpy as np
import pandas as pd
from ete3 import Tree, faces, AttrFace,  TreeStyle

#load tree in newick format
tree = Tree("ML.treefile", format=1)

#defines ancestor of these species as tree outgroup
ancestor = tree.get_common_ancestor("MUSCA_DOMESTICA","GLOSSINA_MORSITANS")
tree.set_outgroup(ancestor)

#load annotated genome info file(file str: specie_sname taxid accession busco_score total_genes)
file = pd.read_csv("./annotatedAssemblies.info.txt", index_col=None, sep='\t', header=None)
#print(file)

#load all genome info file(file str: species_name assembly_name busco_score ref_quality cactus_name) 
genome_info = pd.read_csv("./genome_info.txt", index_col=None, sep='\t', header=None)

#make a list of annotated species
annotated_species = []
for sp in file[0]:
    new_sp = sp.rstrip("\n").replace(" ","_").upper()
    annotated_species.append(new_sp)

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
        else:
            count = count + 1
            leaf.add_features(name=f"Anc{count}")
            
tree_annotation(tree)

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
def get_ref_sp(node, file, busco_relaxation = 0.01):
    #input: node- tree/subtree and file/(dataframe/table) containing species name and busco score
    busco_dict={}
    total_genes_dict={}
    for sp in node.search_nodes(annotated=True):
        sp_name = str(sp).replace("--", "").strip("\n")
        busco_score = file.loc[file[0].str.upper().str.replace(" ","_") == sp_name, 3].item()           
        busco_dict[sp_name] = busco_score
        total_genes = file.loc[file[0].str.upper().str.replace(" ","_") == sp_name, 4].item()
        total_genes_dict[sp_name] = total_genes
    max_busco_sp = max(busco_dict, key=busco_dict.get)
    max_total_genes_sp = max(total_genes_dict, key=total_genes_dict.get)
    new_dict={}
    for sp, busco in busco_dict.items():
        if busco_dict[sp] >= busco_dict[max_busco_sp] - busco_relaxation:
            new_dict[sp] = total_genes_dict[sp]
    return min(new_dict.items(), key=criterion)[0]  #returns a key(ref_sp) having value closest to newValue 
    
#Make a cactus config file containig subtree on line1, Root on line3, Ref_sp on line4, and cactus species name and location
#of their assemblies
def write_treefile_cactus(node, ref_sp, genome_info):
    root_name = node.name
    other_annotated_sp = ""
    for count in range(len(node.search_nodes(annotated=True))):
        new_sps = str(node.search_nodes(annotated=True)[count]).replace("--", "").strip("\n") + ", "
        other_annotated_sp = other_annotated_sp + new_sps
    node.write(format=1, outfile=f"{root_name}_{ref_sp}_tree.nw")
    with open(f"{root_name}_{ref_sp}_tree.nw", 'a') as f:
        print("\n\n#Root: " + str(node.name) +"\n#Reference_species: " + str(ref_sp) \
              +"\n#Annotated_species: " + other_annotated_sp.rstrip(", ") + "\n", file=f)
        for sp in node.get_leaf_names(is_leaf_fn=None):
            assembly = genome_info.loc[genome_info[0].str.upper().str.replace(" ","_").str.replace(r".","_", regex=True) == sp, 1].item()
            cactus_name = genome_info.loc[genome_info[0].str.upper().str.replace(" ","_").str.replace(r".","_", regex=True) == sp, 4].item()
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
                write_treefile_cactus(node, ref_sp, genome_info)                
                #remove already included non_ref species and collapse subtree(node) to ref species
                for sp in node.get_leaf_names(is_leaf_fn=None):
                    if sp != ref_sp:
                        sp = tree.search_nodes(name=sp)[0]
                        sp.delete()
            elif len(node.search_nodes(annotated=True)) > 1:
                #select species with highest busco as a reference
                root_name = node.name
                ref_sp = get_ref_sp(node, file)
                write_treefile_cactus(node, ref_sp, genome_info)
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


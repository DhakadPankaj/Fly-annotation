#!/usr/bin/env python 

#Tree plotting
#Usage: python3 subclades_highlight_plot.py treefile genome_info_file folder_cactus_files

import sys
import os
import random
import colorsys
import itertools
import numpy as np
import pandas as pd
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace

print(sys.argv)      
treefile = sys.argv[1]  
genome_info_file = sys.argv[2]
folder_cactus_files = sys.argv[3]

#load tree in newick format
t_draw = Tree(treefile, format=1)
ancestor = t_draw.get_common_ancestor("MUSCA_DOMESTICA","GLOSSINA_MORSITANS")
t_draw.set_outgroup(ancestor)

#print(file)
genome_info = pd.read_csv(genome_info_file, index_col=None, sep='\t')

#make a list of annotated species
annotated_species = genome_info.loc[genome_info['annotated'] == "Yes", 'hal_name'].tolist()

#make a list of RNA_seq species
rnaseq_species = genome_info.loc[genome_info['RNA_seq'] == "Yes", 'hal_name'].tolist()

#List of species not included in annotation pipeline
outgroup_species = ["GLOSSINA_MORSITANS","MUSCA_DOMESTICA","PAYKULLIA_MACULATA", "POLLENIA_ANGUSTIGENA", "TACHINA_FERA",\
                   "PROTOCALLIPHORA_AZUREA", "STOMORHINA_LUNATA", "SARCOPHAGA_BULLATA", "CIRRULA_HIANS",\
                    "EPHYDRA_GRACILIS","SCATELLA_TENUICOSTA","DIASTATA_REPLETA","BRAULA_COECA"]

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
                leaf.add_features(RNA_seq=True, confidence=1)
        else:
            count = count + 1
            leaf.add_features(name=f"Anc{count}")
            
tree_annotation(t_draw)


def random_color(h=None):
    #Generates a random color in RGB format
    if not h:
        h = random.random()
        s = 0.5
        l = 0.5
    return _hls2hex(h, l, s)

def _hls2hex(h, l, s):
    #converts hls to hex color code
    return '#%02x%02x%02x' %tuple(map(lambda x: int(x*255),
                                        colorsys.hls_to_rgb(h, l, s)))

def random_bgcolor():
    col_list = ["#F8EA94","#F19C92","#92F1C4","#CAEC7E","#BAB2F7","#F7B2DD","#F591A0","#E5FAE3",\
                "#FAE3E3","#C5DEF9","#FAA0FA","#A0CAFA", "#A9971A", "#B36423", "#79B23D", "#46C5BB",\
               "#4684C5", "#7C46C5", "#CC56BE", "#B95A84"]
    return random.choice(col_list)


def layout(node):
    for leaf in t_draw.traverse(): 
        if leaf.name in ref_sp_list:
            #leaf.img_style["fgcolor"] = "#FF0404"
            leaf.img_style["size"] = 13
            #lf.faces_bgcolor = "red"
    if node.is_leaf():
        #if node.name in ref_sp_list and node.name not in rnaseq_species:
         #   N = AttrFace("name", fsize=12, fgcolor = "red")
          #  faces.add_face_to_node(N, node, 0, position="aligned")
        if node.name in rnaseq_species:
            N = AttrFace("name", fsize=12, fgcolor = "blue")
            faces.add_face_to_node(N, node, 0, position="aligned")
        if node.name in alt_rnaseq_sp_list:
            N = AttrFace("name", fsize=12, fgcolor = "red")
            faces.add_face_to_node(N, node, 0, position="aligned")
        if node.name not in rnaseq_species and node.name not in alt_rnaseq_sp_list:
            N = AttrFace("name", fsize=12, fgcolor = "black")
            faces.add_face_to_node(N, node, 0, position="aligned")
            

def get_example_tree(ref_sp):
    # Set dashed blue lines in all leaves
    nst = NodeStyle()
    nst["bgcolor"] = random_bgcolor() 
    #print(ref_sp)
    #print(subtree.get_leaf_names(is_leaf_fn=None))
    n = t_draw.get_common_ancestor(subtree.get_leaf_names(is_leaf_fn=None))
    n.set_style(nst)    
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.scale =  2000
    #ts.mode = "c"
    ts.draw_guiding_lines = True
    ts.legend.add_face(TextFace("Species with their own RNA seq = Blue\nSpecies with RNA seq of closest species= RED\nOther species = BLACK",\
                                fsize=20, fstyle = "italic", fgcolor = "red"), column=1, )
    return t_draw, ts

ref_sp_list=[]
alt_rnaseq_sp_list=[]
for file in os.listdir(folder_cactus_files):
    if file.endswith("_tree.nw"):
        cactus_file = open(f"{folder_cactus_files}/{file}")
        line = cactus_file.readlines()
        subtree = Tree(line[0], format=1)
        ref_sp = line[3].split(" ")[1]
        ref_sp_list.append(ref_sp.rstrip("\n"))
        if len(line[5]) > 1:
            tmp_list = line[5].replace("#Intron_bam: ","").split(",")
            for intron_pair in tmp_list:
                alt_sp = intron_pair.split(";")[0].rstrip("\n").strip(" ")
                alt_rnaseq_sp_list.append(alt_sp)
        t_draw, ts = get_example_tree(ref_sp)
#t_draw.render("subtree_colored.png", w=400, tree_style=ts)
t_draw.show(tree_style=ts)

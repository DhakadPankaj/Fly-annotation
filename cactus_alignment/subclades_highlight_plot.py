#Tree plotting
import os
import random
import colorsys
import itertools
import numpy as np
import pandas as pd
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle


#load tree in newick format
t_draw = Tree("ML.treefile", format=1)
ancestor = t_draw.get_common_ancestor("MUSCA_DOMESTICA","GLOSSINA_MORSITANS")
t_draw.set_outgroup(ancestor)
file = pd.read_csv("./annotatedAssemblies.info.txt", index_col=None, sep='\t', header=None)
#print(file)
genome_info = pd.read_csv("./genome_info.txt", index_col=None, sep='\t', header=None)

annotated_species = []
for sp in file[0]:
    new_sp = sp.rstrip("\n").replace(" ","_").upper()
    annotated_species.append(new_sp)

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
    #print(ref_sp_list)
    for leaf in t_draw.traverse():
        #lf = NodeStyle()
        if leaf.name in ref_sp_list:
            #leaf.img_style["fgcolor"] = "#FF0404"
            leaf.img_style["size"] = 14
            #lf.faces_bgcolor = "red"


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
    ts.show_leaf_name = True
    ts.scale =  2000
    #ts.mode = "c"
    return t_draw, ts

ref_sp_list=[]
for file in os.listdir("./"):
    if file.endswith("_tree.nw"):
        cactus_file = open(file)
        line = cactus_file.readlines()
        subtree = Tree(line[0], format=1)
        ref_sp = line[3].split(" ")[1]
        ref_sp_list.append(ref_sp.rstrip("\n"))
        t_draw, ts = get_example_tree(ref_sp)
#t_draw.render("subtree_colored.png", w=400, tree_style=ts)
t_draw.show(tree_style=ts)


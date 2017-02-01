### Module for Tree Manipulation ###
# -*- coding: latin-1 -*-

# Parse tree
# Input: string tree in input
# Output: tabbed tree

# Use: from pytreetools import *

import sys

def parsed_tree(tree):
    tab_tree = []
    word = ""
    tag = 0
    tag_bracket = 0
    for char in tree:
        if tag == 0:
            if char == "(":
                tag = 1
        if char == "[":
            if len(word) > 0:
                tab_tree.append(word)
            word = char
            tag_bracket = 1
        elif char == "]":
            word = word + char
            tab_tree.append(word)
            word = ""
            tag_bracket = 0
        elif tag_bracket == 1:
            word = word + char
        elif tag == 1:
            if char in ("(", ")", ","):
                tag = 1
                if word != "":
                    tab_tree.append(word)
                word = char
                tab_tree.append(word)
                word = ""
            elif char == ":":
                if word != "":
                    tab_tree.append(word)
                word = char
            elif char == " ":
                tab_tree.append(word)
                word = ""
            elif (char == "." or "_") and (char != ";"):    # Any alphabetical caracter
                word = word + char
            elif char == ";":
                if word != "":
                    tab_tree.append(word)
                #tab_tree.append(char) # Not sure if needed
                break
    return tab_tree


# Definition of the classe Node
class Node:
    """Noeud"""
    def __init__(self):
        self.node = 0         # Number of the node.
        self.parent = 0       # Number of the ancestor of this node
        self.children = []    # List of the children
        self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
        self.duplicate = 0    # Duplicate = 1 if it's a duplicate node
        self.name = ""        # Name of the node: taxa if it's a terminal node, numero if not.
        self.brlen = -1000.0  # Branch length from this node to the ancestor
                              # (-1000 if no branch length)
        self.bootstrap = ""   # Bootstrap (or comment or blabli-blabla)

    def value(self):
        print "Node= "+str(self.node)+" parent= "+str(self.parent)+" Name= "+str(self.name)+ \
              ",\tBranch length = "+str(self.brlen)+",\tBootstrap = "+str(self.bootstrap)+ \
              ",\t  Nodes= "+str(self.children)+", Tip= "+str(self.tip)

    def all_descendant(self): # Find all children from a node
        terminal_nodes  = []
        list_descendant = []
        list_descendant.append(self.node)
        for i in list_descendant:
            if name_object[i].tip == 1:
                terminal_nodes.append(i)
            else:
                for j in name_object[i].children:
                    list_descendant.append(j)
        return terminal_nodes # Return a list
    
    def all_internal_nodes(self): # Find all children from a node
        terminal_nodes = []
        list_nodes = []
        list_nodes.append(self.node)
        for i in list_nodes:
            for j in name_object[i].children:
                list_nodes.append(j)
        return list_nodes # Return a list
    
    def genealogy(self):# Form leaves to root
        ancestres = []
        node = self.node
        while 1:
            if name_object.has_key(node):
                ancestres.append(node)
                node = name_object[node].parent
            else:
                break
        return ancestres # Return a list
    
    def brother(self):    # Find the brother from a node
        for child in name_object[self.parent].children:
            if child == self.node:
                continue
            else:
                brother = child
        return brother

    def split_node(self):
        branches = {}
        children = self.children
        i = 1
        for child in children:
            branches[i]= name_object[child.node].all_descendant()
            i = i+1
        return branches  # Return a dictionnary of lists
    
def read_tree(tree):
    tab_tree = parsed_tree(tree)  # Put the tree in tabbed form 
    node = 0
    level = 0
    actual_node = 0
    terminal_node = 0
    dict_level = {}
    global name_object
    name_object = {}
    for i in range(0, len(tab_tree)):
        item = tab_tree[i]
        #print item
        if item == ";" and inode != 0:  #print "Found ; End procedure";
            break
        elif item == ",":               #print "Found ,"
            continue
        elif item == "(":               #print "Found ("
            node = node+1
            actual_node = node
            if tab_tree[i-1] == "#":
                #print str(nodge)+" duplication!!!!"
                if not name_object.has_key(actual_node):
                    name_object[actual_node] = Node()
                name_object[actual_node].node = actual_node
                name_object[actual_node].duplicate = 1
            level = level+1
            dict_level[level] = node  
        elif item == ")":                #print "Found )"
            actual_node = dict_level[level]
            if not dict_level.has_key(level-1):      # This is the root
                if not name_object.has_key(actual_node):
                    name_object[actual_node] = Node()
                name_object[actual_node].node = actual_node
                name_object[actual_node].parent = -1
                name_object[actual_node].name = "root"
                #if tab_tree[i+1][0] != ":":                              # Need to check
                #    name_object[actual_node].bootstrap = tab_tree[i+1]   # Need to check
                break
            else:
                parent = dict_level[level-1]
                node_fin_parenthese = dict_level[level]
                level = level-1
                if not name_object.has_key(actual_node):
                    name_object[actual_node] = Node()
                name_object[actual_node].node = actual_node
                name_object[actual_node].parent = parent
                name_object[actual_node].name = str(actual_node)
        elif item[0] == ":":  # Branch length found
            if not name_object.has_key(actual_node):
                print "Error: no defined node for branch length"
                break
            name_object[actual_node].brlen = item[1:]
            #print "DISTANCE"
        elif ((item[0] == "#") or (item[0] in str(range(10))) or (item[0] == "[")) \
             and ("_" not in item):  # Bootstrap value found # FOR DRAW NHX
        #elif (((item[0] == "#") or (item[0] in str(range(10)))) or (item[0] == "[")):
            # Bootstrap value found   
        #elif (item[0] == "[") :  # Bootstrap value found in NHX  # FOR EVOLVIEW
            #print "Find duplication node"
            #print "Find bootstrap !!!!!"
            #print item
            if not name_object.has_key(actual_node):           # Need to check
                #print "Error: no defined node for bootstrap   # Need to check
                #break#                                        # Need to check
                name_object[actual_node] = Node()              # Need to check
                name_object[actual_node].node = actual_node    # Need to check
            #print name_object[actual_node].node
            name_object[actual_node].bootstrap = item
            continue
        elif item == "-1":
            continue
        else:    # taxa name
            # Remove bootstrap
            tab_name = item.split("#")
            name = tab_name[0]
            bootstrap = item[len(name):]        
            node = node+1
            actual_node = node
            if not name_object.has_key(node):
                name_object[actual_node] = Node()
            name_object[node].node = actual_node
            name_object[node].name = name
            name_object[node].bootstrap = bootstrap
            name_object[node].parent = dict_level[level]
            name_object[node].tip = 1

    # reached the end of the tree
    
    ### Define children
    for i in name_object:
        list_children=[]
        for l in name_object:
            if name_object[l].parent == name_object[i].node:           
                list_children.append(name_object[l].node)
                #print "Node = "+str(name_object[l].node)+"Parent="+str(name_object[i].node)
            name_object[i].children = list_children
        #print "Child of "+str(name_object[i].node)+" = "+str(name_object[i].children),\
            # name_object[i].tip, name_object[i].name
    return name_object



# Convert taxa name to node
def taxa_to_node(tree, taxa):
    name_object = tree
    node = 0
    for i in name_object:
        if name_object[i].name == taxa:
            node = name_object[i].node
    return node

def node_to_taxa(node):
    taxa = name_object[node].name
    return taxa


# Find Common ancestor node
def common_ancestor(node_list):
    global name_object
    list1 = name_object[node_list[0]].genealogy()
    for node in node_list:
        list2 = name_object[node].genealogy()
        ancestral = []
        for i in list1:
            if i in list2:
                ancestral.append(i)
        list1 = ancestral
    common_ancestor = ancestral[0]
    return common_ancestor # Return a node

# Are all nodes monophyletic ?
def monophyly(node_list):
    global name_object
    is_monophyly = 1
    node_child = ""
    for node in node_list:
        #print name_object[node].genealogy()
        genealogy = name_object[node].genealogy()
        if node_child != "":
            if genealogy[-2] != node_child:
                is_monophyly = 0
        else:
            node_child = genealogy[-2]
        #print name_object[0].brlen
    return is_monophyly # Return 1 (TRUE) or 0 (FALSE)


# Find all internal nodes between two nodes
def connect_nodes(node_list):
    global name_object
    node_root = common_ancestor(node_list)
    list_of_nodes = []
    #list_of_nodes.append(node_root)
    for node in node_list:
        while 1:
            if node == node_root:
                break
            else:
                list_of_nodes.append(node)
                node = name_object[node].parent
    return list_of_nodes

# Get distance from one node to another
def distance_node(node_list):
    global name_object
    if len(node_list) != 2:
        sys.exit()
    node3 = common_ancestor(node_list)
    # Distance from node1 to common ancestor
    distance = 0
    for node in node_list:
        while 1:
            if node == node3:
                break
            else:
                distance = distance+name_object[node].brlen
                node = name_object[node].parent
    return distance

def extract_subtree(node):
    #print "ancetre commun"
    #print node
    # Obtenir la liste de tous les noeuds internes
    #print name_object[node].all_internal_nodes()
    for node1 in name_object[node].all_internal_nodes():
        #print node1
        #print name_object[node1].name
        new_list =  name_object[node].all_internal_nodes()
    new_list.append(1)
    # Supprimer tous les autres
    table = []
    for i in name_object:
        table.append(i)
    for i in table:
        if i not in new_list :
            #print i
            del name_object[i]

    # Move the root.
    #print name_object[1].valeur()
    #print name_object[node].valeur()
    name_object[1].children = name_object[node].children
    for child in name_object[1].children:
        name_object[child].parent = 1
    #print name_object[1].valeur()
    del name_object[node]


# Send a rooted tree and get a unrooted tree
def unroot_tree(tree):
    new_tab_tree = parsed_tree(tree)
    unrooted_tree = []

    level = 1
    tag = 0
    for i in new_tab_tree:
        if i == "(":
            level = level-1
            if ((level != -1) and (tag == 0)) or (tag == 1):
                unrooted_tree.append(i)
        elif i == ")":
            if ((level != -1) and (tag == 0)) or (tag == 1):
                unrooted_tree.append(i)
            if level == -1:
                tag = 1
            level = level+1
        else:
            unrooted_tree.append(i)
    new_tab_tree = "".join(unrooted_tree)
    return new_tab_tree

 
def reroot_tree(new_root):
    ### Rerooting of one tree ###
    ### Get the new and previous distances ###
    old_distance = 0
    old_root = taxa_to_node("root")
    if old_root != 1:
        for i in name_object:
            print name_object[i].valeur()
        sys.exit()

    for i in name_object[old_root].children:
        old_distance = old_distance + name_object[i].brlen
    new_distance = name_object[new_root].brlen
    half_new_distance = name_object[new_root].brlen / 2.0

    ### Get the children form old_root and new_root
    node_right = new_root
    node_left = name_object[new_root].parent
    old_children = name_object[old_root].children
    new_children = [node_left,node_right]

    # Define between the two old root children who's is the parent of who.
    path = name_object[new_root].genealogy()
    new_root_parent = 0
    new_root_child = 0
    new_node2 = 0
    one_of_children = 0
    for node in old_children:
        if node in path:
            new_root_parent = node
            #print "noeud"+str(node)
            new_root_child = name_object[node].brother()
            new_node2 = node
            for node in name_object[node].children:
                if node not in path:
                    one_of_children = node 
            #print "noeud"+str(new_root_child)
    if new_root_parent == 0:
        print "ERROR in the root!!!"
        sys.exit()

    name_object[new_node2].children = [new_root_child,one_of_children]
    name_object[new_root_child].parent = new_root_parent
    name_object[new_root_child].brlen = old_distance
    # Add the distance

    # Inverse the nodes from old_root to new_root
    i=0
    for node in path:
        if node == old_root:
            i = i+1
            continue
        elif node == new_root:
            i = i+1
            continue
        else:
            old_children = name_object[node].children
            new_children123 = []
            for child in old_children:
                if child in path:
                    continue
                else:
                    new_children123.append(child)
            if name_object[node].parent != taxa_to_node("root"):
                new_children123.append(name_object[node].parent)
            if path[i-1] == new_root:
                name_object[node].parent = taxa_to_node("root")
            else:
                name_object[node].parent = path[i-1]
                name_object[node].brlen = name_object[path[i-1]].brlen
            name_object[node].children = new_children123
            #print new_children
            #print name_object[node].parent 
            i=i+1
        #print name_object[node].valeur()

    # Redefine the old_root into new_root
 
    name_object[old_root].children = new_children     # Redefine the children of the root
    name_object[new_root].parent = 1
    name_object[new_root].brlen = half_new_distance   # Add the distance
    name_object[node_left].parent = 1
    name_object[node_left].brlen = half_new_distance  # Add the distance


##def rewrite_tree(node,branch):
##    global new_tree
##    #print list_taxa
##    new_tree = ""
##    new_root = node
##    print node
##    print branch
##    def aff(node):
##        print node
##        #print new_tree
##        #global new_tree
##        print "Noeud="+str(node)
##        print name_object[node].brlen
##        #print name_object[node].valeur()
##        print name_object[node].name
##        distance = str(name_object[node].brlen)
##        distance = distance + ("0"*(7-len(distance)))
##        if node == "":
##            exit
##        if name_object[node].tip == 0:
##            new_tree = new_tree+"("
##            for child in name_object[node].children:
##                aff(child)
##            new_tree = new_tree+")"
##            if node == branch:
##                new_tree = new_tree+"#1"
##            if node != new_root:
##                new_tree = new_tree+":"+distance+","
##        else:
##            new_tree = new_tree+name_object[node].name+":"+distance+","
            
##    aff(node)
##    new_tree = new_tree.replace(",)",")")
##    new_tree = new_tree+";"
##    return new_tree


def write_tree_node():
    global new_tree
    #print list_taxa
    new_tree = ""
    #print node
    #print branch
    node = 1
    branch = 1
    new_root = 1
    def aff(node):
        #print node
        #print new_tree
        global new_tree
        distance = str(name_object[node].brlen)
        distance = str(distance + ("0"*(7-len(distance))))
        if node == "":
            exit
        if name_object[node].tip == 0:
            new_tree = new_tree+"("
            for child in name_object[node].children:
                aff(child)
            new_tree = new_tree+")"
            if node == branch:
                new_tree = new_tree+"#"+str(node)
            if node != new_root:
                new_tree = new_tree+"#"+str(node)+":"+distance+","
        else:
            new_tree = new_tree+name_object[node].name+"#"+str(node)+":"+distance+","

    aff(node)
    new_tree = new_tree.replace(",)",")")
    new_tree = new_tree+";"
    return new_tree

def rewrite_tree(node,branch):
    global new_tree
    new_tree = ""
    new_root = node
    def aff(node):
        global new_tree
        distance = str(name_object[node].brlen)
        distance = str(distance + ("0"*(7-len(distance))))
        if node == "":
            exit
        if name_object[node].tip == 0:
            new_tree = new_tree+"("
            for child in name_object[node].children:
                aff(child)
            new_tree = new_tree+")"
            if node == branch:
                new_tree = new_tree+"#1"
            if node != new_root:
                new_tree = new_tree+name_object[node].bootstrap+":"+distance+","
                #new_tree = new_tree+":"+distance+","
        else:
            new_tree = new_tree+name_object[node].name+name_object[node].bootstrap+":"+distance+","
            #new_tree = new_tree+name_object[node].name+":"+distance+","
    aff(node)
    new_tree = new_tree.replace(",)", ")")
    new_tree = new_tree+";"
    return new_tree

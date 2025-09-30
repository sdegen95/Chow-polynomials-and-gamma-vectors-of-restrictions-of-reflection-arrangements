#########################################################################################################################
# Imports:
#
#     from sage.matroids.advanced import *
#     from collections import OrderedDict
#     from sage.combinat.posets.posets import FinitePoset 
#########################################################################################################################



#########################################################################################################################
# Helper functions for other methods later on
#########################################################################################################################

#-------------------------------------------------------------

def isMinus(s1, s2):
    if len(s1) != len(s2):
        return false
    for a in s1:
        if -a not in s2:
            return false
    return true


#-------------------------------------------------------------

def toString(p):
    if not (isinstance(p, set) or isinstance(p, frozenset)):
        return str(p)
    string  = "{"
    for b in p:
        string = string + str(toString(b))
    return string + "}"


#-------------------------------------------------------------

def merge(partition, b_i, b_j):
    if isMinus(b_i, b_j):
        newpartition = set()
        for b_k in partition:
            if 0 in b_k:
                zero = b_k
                continue
            if b_k != b_i and b_k != b_j:
                newpartition.add(b_k)
        newpartition.add(frozenset(b_i.union(b_j).union(zero)))
    else:
        newpartition = set()
        for b_k in partition:
            if isMinus(b_i, b_k):
                minusb_i = b_k
                continue
            if isMinus(b_j, b_k):
                minusb_j = b_k
                continue
            if b_k != b_i and b_k != b_j:
                newpartition.add(b_k)
        if 0 in b_i or 0 in b_j:
            newblock = b_i.union(b_j.union(minusb_i.union(minusb_j)))
            newpartition.add(frozenset(newblock))
        else:
            newpartition.add(frozenset(b_i.union(b_j)))
            newpartition.add(frozenset(minusb_i.union(minusb_j)))
    return newpartition



#########################################################################################################################
# Methods for computing in intersection lattices of arrangements 
#########################################################################################################################

# Computes the signed set partition lattice, which corresponds to intersection lattice of type $B_n$.
# Input: n = ambient dimension (integer)
# Output: P = intersection lattice of $B_n$ (as FinitePoset object)
#         vertexlabels = Ordered Dictionary of vertex labels, corresponding to the directed graph used for the instantiation of 'P' 

def signed_set_partition_lattice(n): 
    G = DiGraph()
    bottom = set()
    bottom.add(frozenset({0}))
    for i in range(1,n+1):
        bottom.add(frozenset({i}))
        bottom.add(frozenset({-i}))
    G.add_vertex(0)    
    vertexlabels = OrderedDict()
    vertexlabels[0] = bottom
    front = [0]
    nvertices = 0

    while len(front) > 0:
        currentvertex = front.pop()
        partition = vertexlabels[currentvertex]
        for b_i in partition:
            for b_j in partition:
                if b_i == b_j:
                    continue
                

                newpartition = merge(partition, b_i, b_j)
                #newpartition = set()
                #newpartition.add(frozenset(b_i.union(b_j)))
                #for b_k in partition:
                    #if b_k != b_i and b_k != b_j:
                        #newpartition.add(b_k)

                if newpartition not in vertexlabels.values():
                    nvertices += 1
                    vertexnumber = nvertices
                    G.add_vertex(vertexnumber)
                    front.append(vertexnumber)
                    vertexlabels[vertexnumber] = newpartition

                keys = [key for key, val in vertexlabels.items() if val == newpartition]
                newvertex = keys[0]
                G.add_edge(currentvertex, newvertex)
    P = FinitePoset(G)
    
    return P,vertexlabels


# Computes the signed set partition lattices of the intermediate cases $D_ns$, via cutting of some hyperplanes/vertices.
# Input: n = ambient dimension (integer)
#        s = integer between 0 and n
# Output: P = intersection lattice of $D_ns$ (as FinitePoset object)
#         vertexlabels = Ordered Dictionary of vertex labels, corresponding to the directed graph used for the instantiation of 'P'

def d_ns_signed_set_partition_lattice(n,s):
    G = DiGraph()
    bottom = set()
    bottom.add(frozenset({0}))
    for i in range(1,n+1):
        bottom.add(frozenset({i}))
        bottom.add(frozenset({-i}))
    G.add_vertex(0)    
    vertexlabels = OrderedDict()
    vertexlabels[0] = bottom
    front = [0]
    nvertices = 0

    while len(front) > 0:
        currentvertex = front.pop()
        partition = vertexlabels[currentvertex]
        for b_i in partition:
            for b_j in partition:
                if b_i == b_j:
                    continue
                if isMinus(b_i,b_j) and len(b_i)==1 and any([abs(b)>s for b in b_i]):
                    continue
                if 0 in b_i and len(b_i) ==1 and len(b_j)==1 and any([abs(b)>s for b in b_j]):
                    continue
                if 0 in b_j and len(b_j)==1 and len(b_i)==1 and any([abs(b)>s for b in b_i]):
                    continue

                newpartition = merge(partition, b_i, b_j)
                #newpartition = set()
                #newpartition.add(frozenset(b_i.union(b_j)))
                #for b_k in partition:
                    #if b_k != b_i and b_k != b_j:
                        #newpartition.add(b_k)

                if newpartition not in vertexlabels.values():
                    nvertices += 1
                    vertexnumber = nvertices
                    G.add_vertex(vertexnumber)
                    front.append(vertexnumber)
                    vertexlabels[vertexnumber] = newpartition

                keys = [key for key, val in vertexlabels.items() if val == newpartition]
                newvertex = keys[0]
                G.add_edge(currentvertex, newvertex)
    P = FinitePoset(G)
    
    return P,vertexlabels


# Computes the interval atom to top element in a $D_ns$ intersection lattice.
# Input: Po = intersection lattice of a $D_ns$-arrangement (as FinitePoset object)
#        vertexlabels = Ordered Dictionary of vertex labels, corresponding to the directed graph associate to 'Po'
#        k = number of the hyperplane, of whose vertex we want to start with (integer)
# Output: lattice intervall from the atom specify by 'k' to the top element (as FinitePoset object)

def interval_of_lattice(Po,vertexlabels,k):
    n=Po.height()
    bottom = set()
    top=set()
    for i in range(1,n):
        if i != k and i!=-k and i!=0:
            bottom.add(frozenset({i}))
            bottom.add(frozenset({-i}))
            top.add(i)
            top.add(-i)
        else:
            bottom.add(frozenset({-k,0,k}))
            top.add(k)
            top.add(-k)
            top.add(0)
    
    top={frozenset(top)}
    bottom_idx=0
    top_idx=0
    for idx, vertex in vertexlabels.items():
        if vertex == bottom:
            bottom_idx=idx
        if vertex == top:
            top_idx=idx
    return Po.subposet(Po.interval(bottom_idx,top_idx))


# Computes the EL-labeling from Delucchi et.al., for the different signed set partition lattices. 
# We have only sign: zero
# Input: Po = intersection lattice of a $D_ns$-arrangement (as FinitePoset object)
#        VertDict = Ordered Dictionary of vertex labels, corresponding to the directed graph associate to 'Po'
# Output: LabelDict = EL-Labeling dictionary of 'Po'

def EL_Labeling(Po,VertDict):
    CoverRel = [tuple(c) for c in Po.cover_relations()]
    LabelDict = {}
    
    for Rel in CoverRel:
        part_1 = VertDict[Rel[0]]
        part_2 = VertDict[Rel[1]]
        
        
        zero_block_part_1 = [z for z in part_1 if 0 in z][0]
        zero_block_part_2 = [z for z in part_2 if 0 in z][0]
        
        if len(zero_block_part_1) < len(zero_block_part_2): #signed edges
            LabelDict[Rel]=[1,1]
            
        elif len(zero_block_part_1) == len(zero_block_part_2): # unsigned edges
            merging_blocks = list(part_1.difference(part_2))
            merged_blocks = list(part_2.difference(part_1))
            two_merging_blocks = [z for z in merging_blocks if z.issubset(merged_blocks[0])]
            
            rep_1 = min([abs(z) for z in two_merging_blocks[0]])
            sgn_rep_1 = 1 if rep_1 in two_merging_blocks[0] else 0 # 1: contained with positive sign, 0: contained with negative sign
            
            rep_2 = min([abs(z) for z in two_merging_blocks[1]])
            sgn_rep_2 = 1 if rep_2 in two_merging_blocks[1] else 0
            
            if sgn_rep_1 == sgn_rep_2:
                LabelDict[Rel]=[0,max(rep_1,rep_2)]
            elif sgn_rep_1 != sgn_rep_2:
                LabelDict[Rel]=[2,min(rep_1,rep_2)]
                                
    return LabelDict


# Computes all maximal chains in terms of their labelings.
# Input: Po = intersection lattice of a $D_ns$-arrangement (as FinitePoset object)
#        LabelDict = EL-Labeling dictionary of 'Po'
# Output: max_chains_labels = list of all maximal chain, where a maximal chain is list of edge El-Labels 

def EL_Labeling_chains(Po,LabelDict):
    max_chains=Po.maximal_chains()
    max_chains_labels=[]
    for chain in max_chains:
        chain_labels=[]
        for i in range(len(chain)-1):
            chain_labels.append(LabelDict[tuple((chain[i],chain[i+1]))])
        max_chains_labels.append(chain_labels)
    return max_chains_labels


# Counts the number of descents for a single chain.
# (a,b) is a descent if a>=b
# Input: chain  = list of edge El-Labels in the El-labeled intersection lattice of a $D_ns$-arrangement.
# Output: number of descents in the 'chain' (integer) 

def descents_of_chain(chain):  
    des = 0
    for i in range(len(chain)-1):
        if chain[i]>=chain[i+1]:
            des+=1
    return des


# Computes the Chow polynomial for the different signed set partition lattices.
# Input: Po = intersection lattice of a $D_ns$-arrangement (as FinitePoset object)
#        LabelDict = EL-Labeling dictionary of 'Po'
# Output: Chow polynomial 

def Chow_poly(Po,LabelDict):
    n = Po.height()-1
    
    # Compute valid chain for Chow poly
    chain_labels = EL_Labeling_chains(Po,LabelDict)
    valid_chain_labels = []
    for c in chain_labels:
        if c[0]<c[1]:
            temp = True
            for i in range(2,len(c)-1):
                if c[i] >= c[i+1]:
                    if c[i-1] >= c[i]:
                        temp = False
                        break
            if temp == True:
                valid_chain_labels.append(c)

    # Compute the Chow polynomial
    t = var('t')
    Chow_poly = sum([t^(descents_of_chain(c))*(1+t)^(n-1-2*descents_of_chain(c)) for c in valid_chain_labels])
        
    return Chow_poly.expand()


















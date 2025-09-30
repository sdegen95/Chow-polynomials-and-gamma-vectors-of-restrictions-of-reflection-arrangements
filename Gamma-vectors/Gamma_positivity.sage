
#########################################################################################################################
# Imports:
#
# You need the "Create_SimplHyperArr_TopeGraph.sage" file from the GitHub repository:
#
#       https://github.com/Tobias271828/Hamilton-Cycles-in-Simplicial-Arrangements 
#########################################################################################################################



#########################################################################################################################
# Methods for tope graphs, h-vectors and $\gamma$-vectors of a simplicial arrangement
#########################################################################################################################

# Computes all normal vectors of an arrangement.

def normals(arr):
    hyperps = list(arr.hyperplanes())
    normals = [list(h.normal()) for h in hyperps]
    return normals


# Computes the hvector using the topegraph, where the vertices labeled with sign vectors

def hvector(G,d):
    D = []
    for u,v,w in G.edges():
        if u.count(-1) < v.count(-1):
            D.append( (u,v) )
        else:
            D.append( (v,u) )
    D = DiGraph( D) 
    h = [ 0 for i in range(d+1) ]
    for v in D.vertices():
        out_deg = D.out_degree(v)
        h[ out_deg ] += 1
    return h


# Computes the hvector using the topegraph, where the sign vectors for the vertices are stored in a dictionary

def hvector_alt(G,vertInd,d):
    D = []
    for u,v,w in G.edges():
        if len(vertInd[0][u]) < len(vertInd[0][v]):
            D.append( (u,v) )
        else:
            D.append( (v,u) )
    D = DiGraph( D) 
    h = [ 0 for i in range(d+1) ]
    for v in D.vertices():
        out_deg = D.out_degree(v)
        h[ out_deg ] += 1
    return h


# Compute the gamma-vector from the h-vector

def gamma_vector(h, d):
    t = var('t')
    h_polynomial = sum([h[i]* t^i for i in range(d+1)])
        
    gamma = []
    hh = h_polynomial
    while hh != 0:
        hhdeg = hh.degree(t)
        hhleadcoeff = hh.leading_coeff(t)
        gamma.append( hhleadcoeff )
        hh = hh - hhleadcoeff * t^(d-hhdeg) * (1+t)^(d - 2*(d-hhdeg))
        hh = hh.expand()
    print(f"gamma = {gamma}")
    return gamma


# Computes the gamma vector using the h-vector
# solves the problem with 0 in the gamma vector

def gamma_vector_alt(h_vec):
    gam_vec = []
    n = len(h_vec)-1
    d = n//2

    # gamma 0
    gam_vec.append(h_vec[0])

    # gamma 1 to d
    for i in range(1,d+1):
        g_k = h_vec[i]
        for j in range(i):
            g_k -= math.comb(n-2*j,i-j)*gam_vec[j]
        gam_vec.append(g_k)

    return gam_vec



#########################################################################################################################
# Methods fore $D_{n,s}$-arrangements
#########################################################################################################################

# Computes the D_ns arrangement for given d and s where d=n here

def dns_arrangement(d,s):
    var_names = [f'x{i}' for i in range(1,d+1)]
    var_names=tuple(var_names)
    H = HyperplaneArrangements(QQ, names=var_names)
    e = identity_matrix(d).columns()
    PosRoots = [e[i] for i in range(s)]
    PosRoots += [ e[i] - e[j] for i in range(d) for j in range(i+1,d)]
    PosRoots += [ e[i] + e[j] for i in range(d) for j in range(i+1,d)]
    Input = [[tuple(v),0] for v in PosRoots]
    D_ns = H(Input)
    
    return D_ns


# Computes the positive roots of the D_ns arrangement.

def dns_arrangement_posRoots(d,s):
    e = identity_matrix(d).columns()
    PosRoots = [e[i] for i in range(s)]
    PosRoots += [ e[i] - e[j] for i in range(d) for j in range(i+1,d)]
    PosRoots += [ e[i] + e[j] for i in range(d) for j in range(i+1,d)]
    
    return PosRoots


# Computes the topograph for the D_ns arrangements.

def dns_arrangement_topegraph(d,s):
    e = identity_matrix(d).columns()
    PosRoots = [e[i] for i in range(s)]
    PosRoots += [ e[i] - e[j] for i in range(d) for j in range(i+1,d)]
    PosRoots += [ e[i] + e[j] for i in range(d) for j in range(i+1,d)]
    G, vertInd = KonstruiereGraph(PosRoots, d)
    
    return G, vertInd



#########################################################################################################################
# Predefine variables for the $E_7$ and $E_8$ arrangments
#########################################################################################################################


# Global ambient dimension for both cases
d=8


# $E_7$-Arrangement
###################

# E_7 root system
E7_vec_1 = [1,-1,0,0,0,0,0,0]
E7_vec_2 = [1/2,1/2,1/2,1/2,-1/2,-1/2,-1/2,-1/2]
E7_roots = [vector(v) for v in Permutations(E7_vec_1)]
E7_roots += [vector(v) for v in Permutations(E7_vec_2)]


# Hyperplane Arrangement for E_7
var_names = [f'x{i}' for i in range(1,d+1)]
var_names=tuple(var_names)
E7_H = HyperplaneArrangements(QQ, names=var_names)
E7_posRoots=[]
for v in E7_roots:
    if -v not in E7_posRoots:
        E7_posRoots+=[v]     
E7_HA =E7_H([[tuple(v),0] for v in E7_posRoots])


# Restriction of E_7 Arrangement
E7_hyperps = list(E7_HA.hyperplanes())
E7_restriction = E7_HA.restriction(E7_hyperps[0])
E7_restriction



# $E_8$-Arrangement
###################

# E_8 root system
E8_vec_1 = [1,-1,0,0,0,0,0,0]
E8_vec_2 = [-1,-1,0,0,0,0,0,0]
E8_vec_3 = [1,1,0,0,0,0,0,0]
E8_vec_4 = [1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2]          # 0 minus signs
E8_vec_5 = [1/2,1/2,1/2,1/2,1/2,1/2,-1/2,-1/2]        # 2 minus signs
E8_vec_6 = [1/2,1/2,1/2,1/2,-1/2,-1/2,-1/2,-1/2]      # 4 minus signs
E8_vec_7 = [1/2,1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2]    # 6 minus signs
E8_vec_8 = [-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2]  # 8 minus signs
E8_roots = [vector(v) for v in Permutations(E8_vec_1)]
E8_roots += [vector(v) for v in Permutations(E8_vec_2)]
E8_roots += [vector(v) for v in Permutations(E8_vec_3)]
E8_roots += [vector(v) for v in Permutations(E8_vec_4)]
E8_roots += [vector(v) for v in Permutations(E8_vec_5)]
E8_roots += [vector(v) for v in Permutations(E8_vec_6)]
E8_roots += [vector(v) for v in Permutations(E8_vec_7)]
E8_roots += [vector(v) for v in Permutations(E8_vec_8)]


# Hyperplane Arrangement for E_8
var_names = [f'x{i}' for i in range(1,d+1)]
var_names=tuple(var_names)
E8_H = HyperplaneArrangements(QQ, names=var_names)
E8_posRoots=[]
for v in E8_roots:
    if -v not in E8_posRoots:
        E8_posRoots+=[v]
E8_HA =E8_H([[tuple(v),0] for v in E8_posRoots])


# Restriction of E_8 Arrangement to Dimension 7
E8_hyperps = list(E8_HA.hyperplanes())
E8_restriction_dim7 = E8_HA.restriction(E8_hyperps[0])


# Restriction of E_8 Arrangement to Dimension 6
E8_restriction_dim7_hyperps = list(E8_restriction_dim7.hyperplanes())
E8_restriction_dim6_1 = E8_restriction_dim7.restriction(E8_restriction_dim7_hyperps[0])
E8_restriction_dim6_2 = E8_restriction_dim7.restriction(E8_restriction_dim7_hyperps[20])





















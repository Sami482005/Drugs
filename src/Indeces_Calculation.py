import Reading_chem_files as RF
import numpy as np

def getEdgeDensity(mol):
    data = RF.get_molecule_matrices(mol)
    adj = data["adjacency"]
    num_atoms = data["num_atoms"]

    num_edges = 0
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            num_edges += adj[i][j]

    max_edges = num_atoms * (num_atoms - 1) / 2
    return round(num_edges / max_edges, 4) if max_edges else 0

'''
The Wiener Index is defined as the sum of the shortest path distances between all pairs of 
atoms (nodes) in a molecule, where the shortest path between two atoms represents the 
minimum number of bonds (edges) connecting them.

It is calculated by considering the shortest path distances between all pairs 
of atoms (nodes) in the molecule.
Here we will add only half of the matrix, because the distance matrix is symmetric.
'''
def getWeinerIndex(mol):
    data = RF.get_molecule_matrices(mol)
    dist = data["distance"]
    
    N = len(dist)
    wiener_index = 0
    for i in range(N):
        for j in range(i + 1, N):
            wiener_index += dist[i][j]
            #For each pair of nodes (i, j), we retrieve the shortest path 
            #distance dist[i][j] from the distance matrix and add it to wiener_index.
    return wiener_index

'''
Petitjean index is defined the eccentricity of an atom in a
molecular structure as the longest path between that
atom and any other atom in the structure

For that, we are looking to calculate the eccentricities
of each node in the graph. So, we need to calculate
the longest shortest path of every atom.

We define a function LongestShortest_paths(matrix), 
taking matrix as an input

'''

def getPetitjeanIndex(mol):
    data = RF.get_molecule_matrices(mol)
    dist = data["distance"]
    eccentricities = LongestShortest_paths(dist)
    D = max(eccentricities)
    R = min(eccentricities)
    return round((D - R) / R, 4)

'''
Now, in order to find the diameter D and radius R for the calculation of
the petitjean index, we define a function finding_D_R(eccentricities)
diameter D = maximum longest shortest path distance between any 2 atoms in the graph.
radius R = minimum longest shortest path in the graph.
'''
'''
def LongestShortest_paths(matrix):
    N = len(matrix)
    return [max(row[i] for i in range(N) if i != j) for j, row in enumerate(matrix)]
'''
def LongestShortest_paths(matrix):
  
    N = len(matrix)  
    eccentricities = [0] * N   
    for i in range(N):
        Max = 0  
        for j in range(N):
            if i != j:
                Max = max(Max, matrix[i][j])     
        eccentricities[i] = Max
    
    return eccentricities


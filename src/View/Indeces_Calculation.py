import os
import Reading_chem_files as RF

def getEdgeDensity(file_path):
    """
    Calculate the edge density of a molecule from a file.
    
    Parameters:
    - file_path: Path to the input file (sdf or mol).
    
    Returns:
    - edge_density: Calculated edge density.
    """
    # Placeholder for actual edge density calculation
    edge_density = 0.0
    # Implement the logic to read the file and calculate edge density
    return edge_density
def getWeinerIndex(file_path):
    """
    Calculate the Weiner index of a molecule from a file.
    Parameters:
    - file_path: Path to the input file (sdf or mol).
    Returns:
    - weiner_index: Calculated Weiner index.
    """
    # Placeholder for actual Weiner index calculation
    weiner_index = 0.0
    # Implement the logic to read the file and calculate Weiner index
    return weiner_index

def getPetitjeanIndex(file_path):
    '''
    Now, in order to find the diameter D and radius R for the calculation of
    the petitjean index, we define a function finding_D_R(eccentricities)
    diameter D = maximum longest shortest path distance between any 2 atoms in the graph.
    radius R = minimum longest shortest path in the graph.
    '''
    matrix = RF.get_molecule_matrices(file_path)
    # Get the distance matrix from the file
    D = max(LongestShortest_paths(matrix))
    R = min(LongestShortest_paths(matrix))
    petitjean_index = D-R/R
    return petitjean_index


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


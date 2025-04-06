import Reading_chem_files as RF
import numpy as np

def getEdgeDensity(file_path):
    data = RF.get_molecule_matrices(file_path)
    adj = data["adjacency"]
    num_atoms = data["num_atoms"]
    num_edges = np.sum(adj) / 2
    max_edges = num_atoms * (num_atoms - 1) / 2
    return round(num_edges / max_edges, 4) if max_edges else 0

def getWeinerIndex(file_path):
    data = RF.get_molecule_matrices(file_path)
    dist = data["distance"]
    weiner = np.sum(dist) / 2
    return round(weiner, 2)

def getPetitjeanIndex(file_path):
    data = RF.get_molecule_matrices(file_path)
    dist = data["distance"]
    eccentricities = LongestShortest_paths(dist)
    D = max(eccentricities)
    R = min(eccentricities)
    return round((D - R) / R, 4) if R != 0 else "Undefined"

def LongestShortest_paths(matrix):
    N = len(matrix)
    return [max(row[i] for i in range(N) if i != j) for j, row in enumerate(matrix)]

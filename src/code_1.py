import os  # Paths
import networkx as nx  # Graphs
import numpy as np  # Matrices + Arrays
import tkinter as tk
from tkinter import filedialog, messagebox, Label
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from rdkit.Chem import Draw
from PIL import ImageTk

def load_molecules(file_path):
    """
    This function loads molecules from a given file or directory. It can handle:
    - A single .mol file.
    - A single .sdf file containing multiple molecules.
    - A directory containing multiple .mol files.

    """
    mols = []  # Initializing an empty list to store molecules

    # Checking if the provided path is a file
    if os.path.isfile(file_path):
        # If the file is a .mol file
        if file_path.endswith(".mol"):
            mol = read_mol_from_file(file_path)  # Reads the molecule from the file
            if mol:  
                mols.append((os.path.basename(file_path), mol))  # Adds the molecule to the list with its name

        # If the file is a .sdf file (can contain multiple molecules)
        elif file_path.endswith(".sdf"):
            # Reads all molecules from the .sdf file
            for i, mol in enumerate(read_mols_from_sdf(file_path)):
                if mol:  # Ensures that the molecule is valid
                    mols.append((f"Molecule_{i+1}", mol))  

    # If the path is a directory (containing multiple .mol files)
    elif os.path.isdir(file_path):
        # Lists all files in the directory
        for file in os.listdir(file_path):
            # If the file is a .mol file
            if file.endswith(".mol"):
                full_path = os.path.join(file_path, file)  # Gets the full path of the file
                mol = read_mol_from_file(full_path)  
                if mol:
                    mols.append((file, mol))  # Adds the molecule to the list with its filename as the name

    # If the input path is neither a valid file nor a directory, it raises an error
    else:
        raise ValueError("Unsupported input. Provide a .mol, .sdf file, or a directory with .mol files.")

    return mols  # Returns the list of molecules

def mol_to_nx(mol):
    """
    Converts an RDKit molecule into a NetworkX graph. In this graph:
    - Nodes represent atoms.
    - Edges represent bonds between atoms.
    
    Each node stores properties like atomic number and symbol, and each edge stores bond type.

    """
    G = nx.Graph()  # Initialize an empty undirected graph

    # Adds atoms as nodes with properties
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), atomic_num=atom.GetAtomicNum(), symbol=atom.GetSymbol())  # Adds atom properties to nodes

    # Adds bonds as edges with properties
    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=str(bond.GetBondType())  # Adds bond type as an edge property
        )

    return G 


def read_mols_from_sdf(file_path):
    
   # Reads multiple molecules from an .sdf file.
    suppl = SDMolSupplier(file_path)  # Uses RDKit's SDMolSupplier to read molecules from the .sdf file
    return [mol for mol in suppl if mol is not None]  # Returns a list of valid molecules


def read_mol_from_file(file_path):
    #Reads a single molecule from a .mol file using RDKit.

    return Chem.MolFromMolFile(file_path)  # Uses RDKit's Chem.MolFromMolFile to read a single molecule


def get_molecule_matrices(mol):
    """
    Generates and returns matrices representing the adjacency matrix, distance matrix, and number of atoms in the molecule.
    - Adjacency matrix represents the bonds between atoms.
    - Distance matrix represents the shortest path lengths between atoms.

    """
    G = mol_to_nx(mol)  # Converts the molecule to a NetworkX graph
    num_atoms = G.number_of_nodes()  # Gets the number of atoms (nodes) in the graph
    adj_matrix = nx.to_numpy_array(G, dtype=int)  # Converts the graph to an adjacency matrix (a matrix of bonds)

    # Calculates shortest path lengths between all pairs of atoms using NetworkX
    length_dict = dict(nx.all_pairs_shortest_path_length(G))
    dist_matrix = np.zeros((num_atoms, num_atoms), dtype=int)  # Initializes an empty distance matrix

    # Loop through all atom pairs (i, j) and populates the distance matrix with the shortest path lengths.
    for i in range(num_atoms):
        for j in range(num_atoms):
            dist_matrix[i][j] = length_dict[i].get(j, -1)  # If there's no path between atoms i and j, set the distance to -1.

    return {
        "adjacency": adj_matrix,
        "distance": dist_matrix,
        "num_atoms": num_atoms,
        "graph": G
    }


def getEdgeDensity(mol):
    """
    Calculates the edge density of a molecule from the adjacency matrix.

    Edge density is the ratio of the number of edges (bonds) to the maximum possible number of edges in the molecule.
    The calculated edge density is rounded to 4 decimal places.
    """
    data = get_molecule_matrices(mol)  # Gets molecule matrices (adjacency, distance)
    adj = data["adjacency"]
    num_atoms = data["num_atoms"]

    num_edges = np.sum(adj) // 2  # Calculates the number of edges (we divide by 2 to avoid double-counting)
    max_edges = num_atoms * (num_atoms - 1) / 2  # Maximum possible number of edges

    return round(num_edges / max_edges, 4)  # Returns the edge density rounded to 4 decimal places


def getWeinerIndex(mol):
    """
    Calculates the Wiener Index of a molecule from its distance matrix.
    
    The Wiener Index is the sum of the shortest path distances between all pairs of atoms, divided by 2.
    
    """
    data = get_molecule_matrices(mol)  
    dist = data["distance"]

    N = len(dist)  # Number of atoms (nodes)
    wiener_index = 0

    # Sums the shortest path lengths between all pairs of atoms
    for i in range(N):
        for j in range(i + 1, N):  # Avoid double-counting by starting from i + 1
            wiener_index += dist[i][j]

    return wiener_index 


def getPetitjeanIndex(mol):
    """
    Calculating the Petitjean index of a molecule.
    
    The Petitjean index is calculated using the formula (D - R) / R, where D is the maximum eccentricity
    and R is the minimum eccentricity of the atoms in the molecule.
    
    """
    data = get_molecule_matrices(mol)  
    dist = data["distance"]

    # Calculates the eccentricities of the atoms (nodes)
    eccentricities = LongestShortest_paths(dist)

    # Calculating the Petitjean index using the formula (D - R) / R
    D = max(eccentricities)
    R = min(eccentricities)

    return round((D - R) / R, 4)  # Returns the Petitjean index rounded to 4 decimal places


def LongestShortest_paths(matrix):
    """
    Calculate the eccentricity of each node (atom) based on the shortest path matrix.
    
    Eccentricity of a node is defined as the maximum shortest path distance from that node to any other node.

    """
    N = len(matrix) 
    eccentricities = [0] * N  

    # Calculates the maximum shortest path length for each node
    for i in range(N):
        Max = 0
        for j in range(N):
            if i != j:  # Doesn't calculate the path length from a node to itself
                Max = max(Max, matrix[i][j])  # Finds the longest shortest path
        eccentricities[i] = Max  # Stores the eccentricity for node i

    return eccentricities  

class FileUploaderApp(tk.Frame):
    """
    This class represents a Tkinter GUI for uploading files or directories, 
    selecting indices to calculate (Edge Density, Wiener Index, Petitjean Index), 
    and displaying results.
    """
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.parent.title("Indices Calculator")
        self.parent.geometry("500x500")
        self.frames = []
        self.selected_path = None
        
        tk.Label(self, text="Welcome to the Indices Calculator", fg="blue", font=("Times New Roman", 24)).pack(pady=10)
        self.edge_var = tk.BooleanVar()
        self.weiner_var = tk.BooleanVar()
        self.petitjean_var = tk.BooleanVar()


        checkbox_frame = tk.Frame(self)
        checkbox_frame.pack(pady=5)
        tk.Checkbutton(checkbox_frame, text="Edge Density", variable=self.edge_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)
        tk.Checkbutton(checkbox_frame, text="Weiner Index", variable=self.weiner_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)
        tk.Checkbutton(checkbox_frame, text="Petitjean", variable=self.petitjean_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)

        # Buttons for uploading files or directories
        tk.Button(self, text="Upload sdf or mol file", font=("Times New Roman", 18), command=self.upload_file).pack(pady=5)
        tk.Button(self, text="Upload Directory of mol files", font=("Times New Roman", 18), command=self.upload_directory).pack(pady=5)


        self.selected_label = tk.Label(self, text="No file or directory selected", fg="red", font=("Times New Roman", 20))
        self.selected_label.pack(pady=5)

        tk.Button(self, text="Go", font=("Times New Roman", 20), command=self.open_correct_frame).pack(side=tk.RIGHT, padx=10, pady=10)


    def upload_file(self):
        """
        This method handles the file upload action. It opens a file dialog to allow the user to
        select a single file (either .sdf or .mol) from the filesystem.
        
        After selecting the file, it updates the label with the file's name and path.
        """
        file_path = filedialog.askopenfilename(filetypes=[("Chemical Files", "*.sdf *.mol")])  # Opens file dialog
        if file_path:  
            self.selected_path = file_path  # Stores the selected file path
            # Updates the label to show the selected file's name
            self.selected_label.config(text=f"Selected file: {os.path.basename(file_path)}", fg="green", font=("Times New Roman", 18))


    def upload_directory(self):
        """
        This method handles the directory upload action. It opens a directory dialog to allow the user 
        to select a directory containing .mol files.
        
        After selecting the directory, it updates the label with the directory name.
        """
        dir_path = filedialog.askdirectory()  # Opens directory dialog
        if dir_path:  
            self.selected_path = dir_path  # Stores the selected directory path
            # Updates the label to show the selected directory's name
            self.selected_label.config(text=f"Selected directory: {os.path.basename(dir_path)}", fg="green", font=("Times New Roman", 18))


    def open_correct_frame(self):
        """
        This method is triggered when the user clicks the 'Go' button. It checks whether:
        1. At least one index (Edge Density, Wiener Index, Petitjean Index) is selected.
        2. A file or directory has been uploaded.
        
        If both conditions are met, it proceeds to the next frame for results display.
        """
        # Checks if at least one index is selected
        if not (self.edge_var.get() or self.weiner_var.get() or self.petitjean_var.get()):
            messagebox.showerror("Error", "Please select at least one index to calculate.")  # Shows error message if an index wasn't selected
            return

        # Checks if a file or directory is uploaded
        if not self.selected_path:
            messagebox.showerror("Error", "Please upload a file or directory first.") 
            return  

        # Stores the selected indices to be calculated
        self.selected_indices = []
        if self.edge_var.get(): self.selected_indices.append("Edge Density")
        if self.weiner_var.get(): self.selected_indices.append("Weiner Index")
        if self.petitjean_var.get(): self.selected_indices.append("Petitjean Index")

        # Checks if the selected path is valid (either a .mol file, .sdf file, or directory)
        if self.selected_path.endswith(".mol") or self.selected_path.endswith(".sdf") or os.path.isdir(self.selected_path):
            # Creates a new ResultsFrame to show the results
            frame = ResultsFrame(self, self, self.selected_path, self.selected_indices)
            self.show_frame(frame)  


    def show_frame(self, frame):
        """
        This method handles switching between frames in the GUI.
        It destroys the current frame (the fileuploader) and displays the new frame.
        """
        if self.frames:
            old_frame = self.frames.pop()
            old_frame.destroy()
        
        self.frames.append(frame)
        frame.pack(fill="both", expand=True)


class ResultsFrame(tk.Frame):
    """
    This class creates a frame that displays the results after processing the uploaded molecules.
    It shows the calculated indices and their values for each molecule.
    """
    def __init__(self, parent, controller, file_path, selected_indices):
        for widget in parent.winfo_children():
            widget.destroy()

        super().__init__(parent)
        self.controller = controller
        self.file_path = file_path
        self.selected_indices = selected_indices
        self.controller.parent.title("Mol File Results")  # Change window title
        self.controller.parent.geometry("700x800")  # Resize the window

        # Creates a canvas with a scrollbar for displaying results
        canvas = tk.Canvas(self)
        scrollbar = tk.Scrollbar(self, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.pack(fill="both", expand=True)

        molecules = load_molecules(file_path)
        for name, mol in molecules:
            G = mol_to_nx(mol)  # Convert molecule to graph
            # Displays the results for each molecule
            self.display_molecule_results(scrollable_frame, name, G, mol)


    def display_molecule_results(self, frame, name, G, mol):
        """
        Display the results for each molecule: image, name, and calculated indices.
        """
        # Convert RDKit molecule to an image
        img = Draw.MolToImage(mol, size=(600, 600))
        img_tk = ImageTk.PhotoImage(img)

        # Displays image in Tkinter
        img_label = tk.Label(frame, image=img_tk)
        img_label.image = img_tk  # Keep a reference to avoid garbage collection
        img_label.pack(pady=10)

        # Displays the molecule name
        Label(frame, text=f"Molecule: {name}", font=("Times New Roman", 14, "bold")).pack(pady=(5, 10))

        # Calculates selected indices for the molecule
        indices = self.calculate_indices(mol)

        # Displays results for selected indices
        if "Petitjean Index" in self.selected_indices:
            Label(frame, text=f"PetitJean: {indices.get('Petitjean', 'N/A')}").pack()
        if "Weiner Index" in self.selected_indices:
            Label(frame, text=f"Weiner: {indices.get('Weiner', 'N/A')}").pack()
        if "Edge Density" in self.selected_indices:
            Label(frame, text=f"Edge Density: {indices.get('EdgeDensity', 'N/A')}").pack()

        # Adds separator between results
        tk.Label(frame, text="-" * 60).pack(pady=15)


    def calculate_indices(self, mol):
        """
        Calculate the selected indices for a given molecule and return the results.
        """
        results = {}  # Dictionary to store results
        try:
            # Calculates Edge Density if selected
            if "Edge Density" in self.selected_indices:
                results["EdgeDensity"] = getEdgeDensity(mol)
            # Calculate Wiener Index if selected
            if "Weiner Index" in self.selected_indices:
                results["Weiner"] = getWeinerIndex(mol)
            # Calculate Petitjean Index if selected
            if "Petitjean Index" in self.selected_indices:
                results["Petitjean"] = getPetitjeanIndex(mol)
        except Exception as e:
            print(f"Error calculating indices for molecule: {e}")  # Prints error if calculation fails
        return results 


# Create the Tkinter main window and run the application
if __name__ == "__main__":
    root = tk.Tk()  # Create the main Tkinter window
    app = FileUploaderApp(root)  # Pass the root as parent to the FileUploaderApp
    app.pack(fill="both", expand=True)  # Pack the FileUploaderApp to the main window
    root.mainloop()  # Start the Tkinter event loop

import os #paths
import networkx as nx #graphs
import numpy as np #matrices + arrays
import matplotlib.pyplot as plt #drawing
import tkinter as tk
from tkinter import filedialog, messagebox, Label
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def parse_mol_block(mol_lines):
    counts_line = mol_lines[3]
    num_atoms = int(counts_line[0:3])
    num_bonds = int(counts_line[3:6])

    atom_lines = mol_lines[4:4 + num_atoms]
    bond_lines = mol_lines[4 + num_atoms:4 + num_atoms + num_bonds]

    G = nx.Graph()

    for i, line in enumerate(atom_lines):
        x = float(line[0:10])
        y = float(line[10:20])
        z = float(line[20:30])
        element = line[31:34].strip()
        G.add_node(i, element=element, x=x, y=y, z=z)

    for line in bond_lines:
        a1 = int(line[0:3]) - 1
        a2 = int(line[3:6]) - 1
        bond_type = int(line[6:9])
        G.add_edge(a1, a2, bond_type=bond_type)

    # Attempt to extract molecule name
    name = mol_lines[0].strip() if len(mol_lines) > 0 else "Unknown_Molecule"

    return name, G

def parse_mol_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return [parse_mol_block(lines)]

def parse_sdf_file(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    mol_blocks = content.strip().split("$$$$")
    results = []

    for block in mol_blocks:
        lines = [line for line in block.strip().split("\n") if line]
        if lines:
            results.append(parse_mol_block(lines))

    return results

def draw_molecule_image(name, G):
    output_path = "images"
    os.makedirs(output_path, exist_ok=True)
    path = os.path.join(output_path, f"{name}.png")

    pos = nx.spring_layout(G, seed=42)
    plt.figure(figsize=(4, 4))
    nx.draw(G, pos, with_labels=True, node_color='lightgreen', edge_color='gray', node_size=500)
    labels = nx.get_node_attributes(G, 'element')
    nx.draw_networkx_labels(G, pos, labels, font_size=10)
    plt.title(name)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    return path

def get_molecule_matrices(G):
    num_atoms = G.number_of_nodes()
    adj_matrix = nx.to_numpy_array(G, dtype=int)

    length_dict = dict(nx.all_pairs_shortest_path_length(G))
    dist_matrix = np.zeros((num_atoms, num_atoms), dtype=int)

    for i in range(num_atoms):
        for j in range(num_atoms):
            dist_matrix[i][j] = length_dict[i][j] if j in length_dict[i] else -1

    return {
        "adjacency": adj_matrix,
        "distance": dist_matrix,
        "num_atoms": num_atoms,
        "graph": G
    }

def get_bonds_from_adjacency(num_atoms, adjacency_matrix):
    bonds = []
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            if adjacency_matrix[i][j] == 1:
                bonds.append((i, j))
    return bonds

def load_molecules(file_path):
    if file_path.endswith(".mol"):
        return parse_mol_file(file_path)
    elif file_path.endswith(".sdf"):
        return parse_sdf_file(file_path)
    if os.path.isdir(file_path):
        file = parse_into_sdf(file_path)
        return parse_sdf_file(file)
    else:
        raise ValueError("Unsupported file type. Only .mol and .sdf are supported.")
def parse_into_sdf(directory):
    output_file = os.path.join(directory, "combined.sdf")
    mol_files = [f for f in os.listdir(directory) if f.endswith(".mol")]
    with open(output_file, 'w') as outfile:
        for mol_file in mol_files:
            with open(os.path.join(directory, mol_file), 'r') as infile:
                content = infile.read().strip()
                outfile.write(content + "\n$$$$\n")
    return output_file
'''
Edge density was the first index discovered historically.
The edge density of a molecule is defined as the ratio of the 
number of bonds in the molecule to the maximum possible number 
of bonds in a molecule with the same number of atoms.
So, literally that's the formula.
'''
def getEdgeDensity(mol):
    data  = get_molecule_matrices(mol)
    adj = data["adjacency"]
    num_atoms = data["num_atoms"]

    num_edges = 0
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            num_edges += adj[i][j]

    max_edges = num_atoms * (num_atoms - 1) / 2
    return round(num_edges / max_edges, 4)


'''
The Wiener Index is defined as the sum of the shortest path distances between all pairs of 
atoms (nodes) in a molecule, where the shortest path between two atoms represents the 
minimum number of bonds (edges) connecting them.

It is calculated by considering the shortest path distances between all pairs 
of atoms (nodes) in the molecule.
Here we will add only half of the matrix, because the distance matrix is symmetric.
'''
def getWeinerIndex(mol):
    data = get_molecule_matrices(mol)
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
    data = get_molecule_matrices(mol)
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


class FileUploaderApp(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)  # Pass the parent to the superclass
        self.parent = parent
        self.parent.title("Indices Calculator")
        self.parent.geometry("500x500")

        self.frames = []  # Stack to keep track of frames

        # Add a welcome label
        tk.Label(self, text="Welcome to the Indices Calculator", fg="blue", font=("Times New Roman", 24)).pack(pady=10)

        # Allow the user to select which indices they want to calculate
        self.edge_var = tk.BooleanVar()
        self.weiner_var = tk.BooleanVar()
        self.petitjean_var = tk.BooleanVar()

        checkbox_frame = tk.Frame(self)
        checkbox_frame.pack(pady=5)
        tk.Checkbutton(checkbox_frame, text="Edge Density", variable=self.edge_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)
        tk.Checkbutton(checkbox_frame, text="Weiner Index", variable=self.weiner_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)
        tk.Checkbutton(checkbox_frame, text="Petitjean", variable=self.petitjean_var, font=("Times New Roman", 20)).pack(side=tk.LEFT, padx=5)

        tk.Button(self, text="Upload sdf or mol file", font=("Times New Roman", 18), command=self.upload_file).pack(pady=5)
        tk.Button(self, text="Upload Directory of mol files", font=("Times New Roman", 18), command=self.upload_directory).pack(pady=5)

        # By default, always tell the user that no file or directory is selected
        self.selected_label = tk.Label(self, text="No file or directory selected", fg="red", font=("Times New Roman", 20))
        self.selected_label.pack(pady=5)

        tk.Button(self, text="Go", font=("Times New Roman", 20), command=self.open_correct_frame).pack(side=tk.RIGHT, padx=10, pady=10)

    def upload_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Chemical Files", "*.sdf *.mol")])
        if file_path:
            self.selected_path = file_path
            self.selected_label.config(text=f"Selected file: {os.path.basename(file_path)}", fg="green", font=("Times New Roman", 18))


    def upload_directory(self):
        dir_path = filedialog.askdirectory()
        if dir_path:
            self.selected_path = dir_path
            self.selected_label.config(text=f"Selected directory: {os.path.basename(dir_path)}", fg="green", font=("Times New Roman", 18))

    def open_correct_frame(self):
        # Error Handling
        if not (self.edge_var.get() or self.weiner_var.get() or self.petitjean_var.get()):
            messagebox.showerror("Error", "Please select at least one index to calculate.")
            return

        if not self.selected_path:
            messagebox.showerror("Error", "Please upload a file or directory first.")
            return

        # Load the indices selected by the user
        self.selected_indices = []
        if self.edge_var.get(): self.selected_indices.append("Edge Density")
        if self.weiner_var.get(): self.selected_indices.append("Weiner Index")
        if self.petitjean_var.get(): self.selected_indices.append("Petitjean Index")
        
        # There are two frames we prepared. one for only 1 molecule (mol file) and one for multiple molecules
        # (sdf file or directory of mol files)
        if self.selected_path.endswith(".mol") or self.selected_path.endswith(".sdf") or os.path.isdir(self.selected_path):
            frame = ResultsFrame(self, self, self.selected_path, self.selected_indices)
            self.show_frame(frame)

    def show_frame(self, frame):
        # Hide the current frame
        if self.frames:
            old_frame = self.frames.pop()
            old_frame.destroy()  # Destroy the old frame
        # Add the new frame to the stack and display it
        self.frames.append(frame)
        frame.pack(fill="both", expand=True)

class ResultsFrame(tk.Frame):
    def __init__(self, parent, controller, file_path, selected_indices):
        # Clear all previous widgets
        for widget in parent.winfo_children():
            widget.destroy()

        super().__init__(parent)
        self.controller = controller
        self.file_path = file_path  
        self.selected_indices = selected_indices  

        self.controller.parent.title("Mol File Results")
        self.controller.parent.geometry("700x800")

        # Set up scrollable canvas
        canvas = tk.Canvas(self)
        scrollbar = tk.Scrollbar(self, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.pack(fill="both", expand=True)

        # Load molecules and display
        molecules = load_molecules(file_path)
        for name, G in molecules:
            self.display_molecule_results(scrollable_frame, name, G)

    def display_molecule_results(self, frame, name, G):
        # --- Draw molecule ---
        fig, ax = plt.subplots(figsize=(5, 5))
        pos = nx.spring_layout(G, seed=42)
        nx.draw(G, pos, ax=ax, with_labels=True, node_color='skyblue', edge_color='gray', node_size=700)
        labels = nx.get_node_attributes(G, 'element')
        nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=10)
        ax.set_title(name)
        ax.axis('off')

        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(pady=10)
        canvas.draw()

        # --- Molecule name ---
        Label(frame, text=f"Molecule: {name}", font=("Times New Roman", 14, "bold")).pack(pady=(5, 10))

        # --- Compute indices ---
        indices = self.calculate_indices(G)

        # --- Show indices ---
        if "Petitjean Index" in self.selected_indices:
            Label(frame, text=f"PetitJean: {indices.get('Petitjean', 'N/A')}").pack()
        if "Weiner Index" in self.selected_indices:
            Label(frame, text=f"Weiner: {indices.get('Weiner', 'N/A')}").pack()
        if "Edge Density" in self.selected_indices:
            Label(frame, text=f"Edge Density: {indices.get('EdgeDensity', 'N/A')}").pack()

        tk.Label(frame, text="-"*60).pack(pady=15)  # Divider

    def calculate_indices(self, G):
        results = {}
        if "Edge Density" in self.selected_indices:
            results["EdgeDensity"] = getEdgeDensity(G)
        if "Weiner Index" in self.selected_indices:
            results["Weiner"] = getWeinerIndex(G)
        if "Petitjean Index" in self.selected_indices:
            results["Petitjean"] = getPetitjeanIndex(G)
        return results
    

if __name__ == "__main__":
    root = tk.Tk()  # Create the main Tkinter window
    app = FileUploaderApp(root)  # Pass the root as parent
    app.pack(fill="both", expand=True)
    root.mainloop()  # Start the Tkinter event loop


# 🌟 Welcome to Our Drug Design Project <img src="assets/medicon.jpg" width="40"/>


## 🧪 Objectives

- [x] Create a menu page where users can upload `.mol`, `.sdf` files, or directories of `.mol` files  
- [X] Calculate the Topological Indices:  
  Users can choose between the **Petitjean Index**, **Wiener Index**, and **Edge Density**
- [X] Display the 2D structure of the molecule with its corresponding index values  

---

## 🧰 Requirements

- Python 3  
- NetworkX (for Graph Theory)  
- Tkinter (for GUI)
- PIL / image viewer (for enhanced visuals)

---

## 📁 File System

### `Src`  
**Contains all backend logic and data processing functions**  
- One python file [Drug_project.py](src/Drug_project.py) that contains both the fileuploader class and the results class, in addition to all the needed methods to complete this project.
---

### `Test_files/`  
**Includes sample `.sdf` and `.mol` files**  
- Used to test uploading functionality, directory handling, and the robustness of the index calculators.
---


## 🚀 Getting Started

```bash
    git clone https://github.com/Sami482005/drug-design-project.git
    cd drug-design-project
    python Drug_project.py
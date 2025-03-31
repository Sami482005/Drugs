# Welcome to Our Drug Design Project <img src='assets/medicon.jpg' width=40>
## Objectives:
- [X] Create a menu page where users can upload mol files, sdf files, or directories of mol files
- [ ] Calculate the Topological Indeces.
Users get to choose between the **Petitjean Index**, **Weiner Index**, and **Edge Density**
- [ ] The Final Menu will produce the 3D structure of the molecule described by the file, and print next to it the value of the Indeces

## Requirements
Python3
Graph Theory
Indeces Calculations

## File_System
### Src/Control:
***Contains all the functions that happen in the backend of the program***
[Reading_chem-files.py](src/Control/Reading_chem_files.py) allows the program to read the mol or sdf file to transform it into an image on the final page.
[Indeces_Calculations.py](src/Control/Indeces_Calculation.py) contains the functions that will calculate the indeces of the molecules.

### Src/View:
***Contains the GUI menus for the program***
#### [First_Page.py](src/View/First_Page.py) 
Shows the main page that will appear to the user once they activate the program. 
Must contain a button "Done" that will take the user to the next page that will contain the results.
#### Results Page:
If only one mol file was selected then the structure would appear in the center of the page, with the name underneath the structure image. With the three values of the indeces displayed beneath it.
If one sdf file or a directory of mol files was selected then the result file would display one molecule's result on the result page. There will be a next and previous page that will display the results of the previous or next molecule in the sdf file or the directory.

### Test_files
Contains a combination of sdf and mol files. Will also test uploading the directory itself to the program to see if it will classify the mol files and if the directory branch of the program works properly.
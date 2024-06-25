# Batoid_Fins
Script and files associated with manuscript on Batoids' pectoral fin shape evolution

Script SelBat_280524.R Script Batoids_310524.R contains all the functions and notes produce the results in the manuscript and supplementary figures/tables.

GM_BatSel folder contains the landmark files, curves and trees files to use with the script files.

PhyloBat  Folder:
The TNT file of the matrix along with the commands used for the analysis. To run the analysis, add the location of the TNT matrix file next to the “cdir” command. Use the command “procedure” in the TNT terminal to run the analysis with the command file (Bat_phylo_analysis_TNTcommands.tnt and Matrix.tnt).
The trees estimated in the Parsimony analysis (Bat_TNT.tre).
The matrix and tree files in NEXUS format needed to run the PAUP analysis along with the commands needed to perform the branch length estimation. Same as in the TNT command file add the location of the NEXUS matrix file next to the “cd” command. Use the command “execute” in the PAUP terminal to run the analysis with the command file (Bat_paup_matrix_trees.nex and Paup_commands.paup)
The log file for both TNT and PAUP analyses (TNT.log and PAUP.log).

# PDBdata
Extract necessary data from PDB_cif
#
cif2pdb.py is the main module that converts PDB mmCIF format file to old PDB format file, when it is used in stand alone. The function preprocess() checks errors in mmCIF file and the function readCIF converts the input file to pdb dictionary. The functions writeHEADER(), writeATOM(), writeHETATM(), and writeEND() make an output in old PDB format.
#
Old but still in use software only accept PDB format file and does not accept mmCIF file. In that case, this converter may solve the problem.
#
The good point of the program is that the program does not need any other extra library except the default ones.

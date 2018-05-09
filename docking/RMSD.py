import numpy as np
import parser_pdb as pdb
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]


# calcule le rmsd a partir d'une liste de distance entre atomes
def rmsd(distances):
    return np.sqrt((np.array(distances) ** 2).mean())

# prends 2 liste des coordonnes des atomes CA concatenee a la suite (x,y,z,x,y,z,...)
# renvoi une liste des distances
def distances(CA1, CA2):
	distanceList = []
	for i in range(0,len(CA1),3):
		distanceList.append(np.sqrt((CA2[i]-CA1[i])**2+(CA2[i+1]-CA1[i+1])**2+(CA2[i+2]-CA1[i+2])**2))
	return distanceList



######## Ex d'Utilisation

#PDB1 = pdb.parse_pdb(file1)
#PDB2 = pdb.parse_pdb(file2)
#coordCa1 = PDB1.getCA()
#coordCa2 = PDB2.getCA()
#distances = distances(coordCa1, coordCa2)
#print(rmsd(distances))





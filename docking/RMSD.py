import parser_pdb as pdb
import sys



# calcule le rmsd a partir d'une liste de distance entre atomes
def rmsd(distances):
    n = len(distances)
    return (sum(map(lambda x: x ** 2, distances)) / n) ** 0.5

# prends 2 liste des coordonnes des atomes CA concatenee a la suite (x,y,z,x,y,z,...)
# renvoi une liste des distances
def distances(CA1, CA2):
	distanceList = []
	for i in range(0,len(CA1),3):
		distanceList.append(((CA2[i]-CA1[i])**2+(CA2[i+1]-CA1[i+1])**2+(CA2[i+2]-CA1[i+2])**2) ** 0.5)
	return distanceList



def rmsdLigant(file_Ligand_Natif, file_Ligand_Soluce):

	PDB1 = pdb.parse_pdb(file_Ligand_Natif)
	PDB2 = pdb.parse_pdb(file_Ligand_Soluce)
	coordCa1 = PDB1.get_CA()
	coordCa2 = PDB2.get_CA()
	distance = distances(coordCa1, coordCa2)
	return rmsd(distance)

def rmsdComplexe(file_Ligand_Natif, file_Ligand_Soluce, file_Rec_Natif):

	PDB1 = pdb.parse_pdb(file_Ligand_Natif)
	PDB2 = pdb.parse_pdb(file_Ligand_Soluce)
	PDB_REC = pdb.parse_pdb(file_Rec_Natif)
	coordCa_Ligand_Natif = PDB1.get_CA()
	coordCa_Ligand_Soluce = PDB2.get_CA()
	coordCa_REC = PDB_REC.get_CA()
	complexe = coordCa_REC + coordCa_Ligand_Natif
	solution = coordCa_REC + coordCa_Ligand_Soluce
	distance = distances(complexe, solution)
	return rmsd(distance)

def rmsdInterface(file_Ligand_Natif, file_Ligand_Soluce, listeIndices):

	PDB1 = pdb.parse_pdb(file_Ligand_Natif)
	PDB2 = pdb.parse_pdb(file_Ligand_Soluce)
	coord1 = PDB1.get_CA_Elems(listeIndices)
	coord2 = PDB1.get_CA_Elems(listeIndices)
	distanceMatrix = distances(coord2, coord2)
	return rmsd(distanceMatrix)



# liste = [0,56] affichera le premier et le 57e element 









from Protein import Pdb
from constants import dcharge, dvdw, depsilon
from parser_pdb import parse_pdb


"""
According to Cronell, 95,
Aij =  epsilon_ij* R_ij*^8
B_ij = 2epsilon_ij* . R_ij*^6
R_ij* = Ri* + Rj*
with R* VdW Radius
epsilon_ij* = (epsilon_i epsilon_j)^(1/2)
"""

class Cornell_calc(object):
    """
    Handles the computation of the free terms in Cornell equations.

    It has an atribute `pdb` wich is a :python:class:`Pdb` object. 
    """    
    const_f = 332.0522    

    def __init__(self, pdb):
        self.pdf = pdb
        self.atom = pdb.get_atom_list()
        self.epsilon = {}
        self.q = {}
        self.r = {}
 
    def get_Eij(self, i, j):
        """Return the free term Eij in the Cornell equation."""

        Aij = self.get_Aij(i, j)
        Bij = self.get_Bij(i, j)
        dist = self.get_dist(i, j)
        qi = self.get_q(i)
        qj = self.get_q(j)
        return Aij / (dist ** 8) - Bij / (dist ** 6) + self.const_f * qi * qj / (20 * dist)     
 
    def get_Aij(self, i, j):
        """Return the free term Aij in the Cornell equation."""
        epsilon_ij = (self.get_epsilon(i) * self.get_epsilon(j))**0.5
        r_ij  = self.get_r(i) + self.get_r(j)
        return epsilon_ij * (r_ij ** 8)     

    def get_Bij(self, i, j):
        """Return the free term Bij in the Cornell equation."""
        epsilon_ij = (self.get_epsilon(i) * self.get_epsilon(j))**0.5
        r_ij  = self.get_r(i) + self.get_r(j)
        return 2 * epsilon_ij * (r_ij ** 6)     

    def get_dist(self, i, j):
        """Return the distances between atom i and atom j."""
        x1, y1, z1 = self.atom[i].getCoord()
        x2, y2, z2 = self.atom[j].getCoord()
        return ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)**0.5     

    def get_epsilon(self, i):
        """Return the espsilon term for atom i."""
        if i not in self.epsilon:
            self.epsilon[i] = None
            self.epsilon[i] = self.get_epsilon(i)
        elif self.epsilon[i] != None:
            return self.epsilon[i]
        
        atomi = self.atom[i].get_context()["atom_type"]
        res_type = self.atom[i].get_context()["residue"]
        if self.is_nter(i):
            if atomi == "H1":
                return 0.0157     
            elif atomi == "H2":
                return 0.0157     
            elif atomi == "H3":
                return 0.0157     
            elif atomi == "N":
                return  0.17     
            elif atomi == "CA":
                return 0.1094     
            elif atomi == "HA":
                return 0.0157     
            else:
                return  depsilon[res_type][atomi]     
        elif self.is_cter(i):
            if atomi == "CA" :
                return 0.1094     
            elif atomi == "C" :
                return 0.0860     
            elif atomi == "O" :
                return 0.2100     
            elif atomi == "OXT" :
                return 0.2100     
            else:
                return  depsilon[res_type][atomi]

        else:
            return  depsilon[res_type][atomi]
              
    def get_q(self, i):
        """Return the charge of atom i."""
        if i not in self.q:
            self.q[i] = None
            self.q[i] = self.get_q(i)
        elif self.q[i] != None:
            return self.q[i]

        atomi = self.atom[i].get_context()["atom_type"]
        res_type = self.atom[i].get_context()["residue"]
        if self.is_nter(i):
            if atomi == "H1" :
                return 0.1984     
            elif atomi == "H2" :
                return 0.1984     
            elif atomi == "H3" :
                return 0.1984        
            elif atomi == "N" :
                return 0.1592     
            elif atomi == "CA" :
                return 0.0221     
            elif atomi == "HA" :
                return 0.116     
            else:
                return dcharge[res_type][atomi]
        elif self.is_cter(i):
            if atomi == "CA" :
                return -0.2493     
            elif atomi == "C" :
                return 0.7231     
            elif atomi == "O" :
                return -0.7855     
            elif atomi == "OXT" :
                return -0.7855     
            else:
                return dcharge[res_type][atomi]
        else:
            return dcharge[res_type][atomi]
                     

    def get_r(self, i):
        """Return the Van der Wals radius of atom i."""
        if i not in self.r:
            self.r[i] = None
            self.r[i] = self.get_r(i)
        elif self.r[i] != None:
            return self.r[i]

        atomi = self.atom[i].get_context()["atom_type"]
        res_type = self.atom[i].get_context()["residue"]    
        if self.is_nter(i):
            if atomi == "H1" :
                return 0.6     
            elif atomi == "H2" :
                return 0.6     
            elif atomi == "H3" :
                return 0.6     
            elif atomi == "N" :
                return 1.875     
            elif atomi == "CA" :
                return 1.9080     
            elif atomi == "HA" :
                return 1.1     
            else:
                return dvdw[res_type][atomi]
        elif self.is_cter(i):
            if atomi == "CA" :
                return 1.9080     
            elif atomi == "C" :
                return 1.9080     
            elif atomi == "O" :
                return 1.6612     
            elif atomi == "OXT" :
                return 1.6612     
            else:
                return dvdw[res_type][atomi]
        else:
            return dvdw[res_type][atomi]
                 
    def is_nter(self, i):
        """Return True if atom i belongs to the N-ter residue."""
        residue = self.atom[i].get_residue()
        if residue._resNum != 1:
            return False    
        residue = self.atom[i].get_residue()
        return ("H1" in residue.keys() or "H2" in residue.keys() or "H3" in residue.keys())    

    def is_cter(self, i):
        """Return True if atom i belongs to the C-ter residue."""
        residue = self.atom[i].get_residue()
        if residue._resNum == 1:
            return False    
        residue = self.atom[i].get_residue()
        return "OXT" in residue.keys()    


def comp_score(receptor_file, ligand_file):
    """Compute the Cornell score for a docking configuration."""

    # Parse the two files
    pdb_receptor = parse_pdb(receptor_file)
    pdb_ligand = parse_pdb(ligand_file)

    # Create a pdb object for the docking
    pdb_docking = Pdb()
    for chainID in pdb_receptor.keys():
        pdb_docking.add_chain("R" + chainID, pdb_receptor[chainID])
    for chainID in pdb_ligand.keys():
        pdb_docking.add_chain("L" + chainID, pdb_ligand[chainID])

    # Create the Cornell ojbect for the docking
    cornell_dock = Cornell_calc(pdb_docking)

    # Split the list of atom's index between ligand and receptor
    last_receptor = len(pdb_receptor.get_atom_list())
    last_ligand = len(pdb_ligand.get_atom_list())
    receptor_idx = range(last_receptor)
    ligand_idx = range(last_receptor, last_receptor + last_ligand)

    #compute the score
    score = 0
    for i in receptor_idx:
        for j in ligand_idx:
            score += cornell_dock.get_Eij(i, j)

    return score




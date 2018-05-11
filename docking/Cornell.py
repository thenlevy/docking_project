from Protein import Pdb
from constants import dcharge, dvdw, depsilon
from parser_pdb import parse_pdb


"""
According to Cronell, 95,
Aij =  epsilon_ij* R_ij*^12
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
 
    def get_Eij(self, i, j):
        """Return the free term Eij in the Cornell equation."""

        Aij = self.get_Aij(i, j)
        Bij = self.get_Bij(i, j)
        dist = self.get_dist(i, j)
        qi = self.get_q(i)
        qj = self.get_q(j)
        return Aij / (dist ** 8) - Bij / (dist ** 6) + f * qi * qj / (20 * dist)     
 
    def get_Aij(self, i, j):
        """Return the free term Aij in the Cornell equation."""
        epsilon_ij = (self.get_epsilon(i) * self.get_epsilon(j))**0.5
        r_ij  = self.get_r(i) + self.get_r(j)
        return epsilon_ij * (r_ij ** 12)     

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
            atomi = atom[i].get_context()["atom_type"]
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

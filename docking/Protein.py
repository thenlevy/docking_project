from numpy import sqrt
class Pdb(dict):
    """
    Represent a complete protein.
    
    This class can be used as a dictionary, that maps chain identifier to
    chain objects.
    """
    def __init__(self, name="pdb"):
        super().__init__()
        self._name = name
        self._chains = []


    def __setitem__(self, key, value):
        if key in self._chains:
            warn_mess = ("Warning, adding " + key + " which is already in "
                         + self._name)
            print(warn_mess, sys.stderr)
        else:
            self._chains.append(key)
        super().__setitem__(key, value)


    def add_chain(self, chain_name, chain):
        """Add a chain object to the protein."""
        self.__setitem__(chain_name, chain)


    def keys(self):
        """Return the list of chain identifiers."""
        return self._chains

    
    def values(self):
        """Return the list of chain objects."""
        return [self[k] for k in self._chains]



    def write(self):
        """Call successively the `write` method of all chain objects."""
        for chain in self.values():
            chain.write()

    def getCA(self):
        """Return the list of atom objects corresponding to the alpha carbon of
        all the residues."""
        out = []
        for chain in self.values():
            out += chain.getCA()
        return out


    def get_atom_list(self):
        """Return a list of atom objects representing all the atoms of the
        protein."""
        ret = []
        for c in self.values():
            ret += c.get_atom_list()
        return ret 



class Chain(dict):
    """
    Represent a chain of a protein.

    This class can be used as a dictionary that maps residue identifier to 
    Residue object.
    """
    def __init__(self, name=""):
        super().__init__()
        self._residue_list = []
        self._name = name
        self._nb_atom = 0
        self._pdb = None

    def __setitem__(self, key, value):
        if key in self._residue_list:
            warn_mess = ("Warning, adding " + key + " which is already in "
                         + self._name)
            print(warn_mess, sys.stderr)
        else:
            self._residue_list.append(key)
        super().__setitem__(key, value)


    def add_residue(self, res_name, res):
        """Add a residue object to the chain."""
        self.__setitem__(res_name, res)
        res.set_chain(self)
        res.set_residue_num()
    

    def keys(self):
        """Return the list of residue identifiers."""
        return self._residue_list


    def values(self):
        """Return the list of residue objects."""
        return [self[k] for k in self._residue_list]


    def set_pdb(self, pdb):
        """Specify the Pdb object to which the chain belongs."""
        self._pdb = pdb


    def provide_residue_num(self):
        """Provide a unique number to a residue object."""
        self._nb_atom += 1
        return self._nb_atom

    def write(self):
        """Call successively the `write` method of all Residue objects."""
        for res in self.values():
            res.write()

    def getCA(self):
        """Return the list of atom objects reprensenting the alpha carbon of all
        the chain's residues"""
        out = []
        for res in self.values():
            out += res.getCA()
        return out

    def get_chain_ID(self):
        """Return the identifier of the chain."""
        return self._name


    def get_atom_list(self):
        """Return the list of the atom objects representing all the atom
        contained in the chain."""
        ret = []
        for r in self.values():
            ret += r.values()
        return ret
            

class Residue(dict):
    """Represent a residue.

    This class can be used as a dictionary maping atom identifiers to atom objects.
    """
    def __init__(self, name=""):
        super().__init__()
        self._atom_list = []
        self._name = name
        self._chain = None
        self._nb_atom = 0


    def __setitem__(self, key, value):
        if key in self._atom_list:
            warn_mess = ("Warning, adding " + key + " which is already in "
                         + self._name)
            print(warn_mess, sys.stderr)
        else:
            self._atom_list.append(key)
            value.set_residue_context(self._name, key)
        super().__setitem__(key, value)


    def add_atom(self, atomtype, atom):
        """Add an atom object."""
        self.__setitem__(atomtype, atom)
        atom.set_residue(self)


    def keys(self):
        """Return the list of atom identifier."""
        return self._atom_list


    def values(self):
        """Return the list of atom objects."""
        return [self[k] for k in self._atom_list]


    def set_chain(self, chain):
        """Specify the chain to which to residue belongs."""
        self._chain = chain


    def set_residue_num(self):
        """Get the number of the residue in its chain."""
        self._resNum = self._chain.provide_residue_num()
        

    def write(self):
        """Successively call the `write` method of all atom objects."""
        for atomtype in self.keys():
            atom = self[atomtype]
            atom.write(atomtype, self._name, self._resNum,
                       self._chain.get_chain_ID())

    def getCA(self):
        """Return the atom object representing the alpha carbon of the residue."""
        out = []
        atom = self['CA']
        out+=atom.getCoord()
        return out

    def get_mass_center(self):
        """Return the coordonates of the center of mass of the residue.

        The weight of the atoms are not taken into consideration.
        """
        x_center = 0
        y_center = 0
        z_center = 0
        for atom in self.values():
            x, y, z = atom.getCoord()
            x_center += x
            y_center += y
            z_center += z
        nb_atom = len(self.values())
        x_center /= nb_atom
        y_center /= nb_atom
        z_center /= nb_atom
        return (x_center, y_center, z_center)

        
class Atom(object):
    """Represent an atom."""
    def __init__(self, x, y, z, identifier, symbol):
        self.x = x
        self.y = y
        self.z = z
        self._identifier = identifier
        self._residue = None
        self._symbol = symbol
        self._context = None

    
    def set_residue(self, residue):
        """Specify the chain to which the atom belongs."""
        self._residue = residue

         
    def write(self, atomtype, resName, resSeq, chain_ID):
        """Print the pdb reprensentation of the atom."""
        out = "ATOM  "
        out += (" " * 5 + str(self._identifier))[-5:] # Atom serial number
        out += "  "
        out += (self._context["atom_type"] + " " * 4)[:4] # Atom name
        out += " " # Alternate location indicator
        out += (self._context["residue"] + " " * 3)[:3] # Residue name
        out += "  "
        out += chain_ID # chain identifier
        out += (" " * 4 + str(resSeq))[-4:] # Residue sequence number
        out += "   "
        out += " " # Code for insertion of residues
        out += (" " * 8 + "{0:.3f}".format(self.x))[-8:] # X coordiates
        out += (" " * 8 + "{0:.3f}".format(self.y))[-8:] # Y coordiates
        out += (" " * 8 + "{0:.3f}".format(self.z))[-8:] # Z coordiates
        out += " " * 6 # Occupancy
        out += " " * 6 # Temperature facotr
        out += " " * 8
        out += (" " * 2 + self._symbol)[-2:] # Element symbol (right justified)
        out += " " # Charge
        print(out)



    def set_residue_context(self, res_name, atom_type):
        """Provide information about the residue to which the atom belongs."""
        self._context = {"residue": res_name, "atom_type": atom_type}


    def get_context(self):
        """Return information about the residue to which the atom belongs."""
        if self._context["residue"] == "HIS":
            if  "HD1" in self._residue.keys():
                self._context["residue"] = "HID"
            else:
                self._context["residue"] = "HIE"
        return self._context

    
    def getCoord(self):
        """Return the coordonates of the atom."""
        out = []
        out.append(self.x) # X coordiates
        out.append(self.y) # Y coordiates
        out.append(self.z) # Z coordiates
        return out


    def get_residue(self):
        """Return the Residue to which the Atom belongs."""
        return self._residue


def distance_residues(r1, r2):
    """Return the distance between two residues."""

    center_1 = r1.get_mass_center()
    center_2 = r2.get_mass_center()
    return sqrt(sum(map(lambda t: (t[0] - t[1])**2, zip(center_1, center_2))))


def interface(chain1, chain2):
    """Return the list of keys of residuals at the interface of chain1 and
    chain2."""
    dist_max = 6
    idx_chain1 = set()
    idx_chain2 = set()
    for k1 in chain1.keys():
        for k2 in chain2.keys():
            if distance_residues(chain1[k1], chain2[k2]) < dist_max:
                idx_chain1.add(k1)
                idx_chain2.add(k2)
    return (idx_chain1, idx_chain2)

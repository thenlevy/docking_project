class Pdb(dict):
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
        self.__setitem__(chain_name, chain)


    def keys(self):
        return self._chains

    
    def values(self):
        return [self[k] for k in self._chains]


    def write(self):
        for chain in self.values():
            chain.write()

    def getCA(self):
        out = []
        for chain in self.values():
            out += chain.getCA()
        return out


    def get_atom_list(self):
        ret = []
        for c in chain:
            ret += c.get_atom_list()
        
    


class Chain(dict):
    def __init__(self, name=""):
        super().__init__()
        self._residual_list = []
        self._name = name
        self._nb_atom = 0
        self._pdb = None

    def __setitem__(self, key, value):
        if key in self._residual_list:
            warn_mess = ("Warning, adding " + key + " which is already in "
                         + self._name)
            print(warn_mess, sys.stderr)
        else:
            self._residual_list.append(key)
        super().__setitem__(key, value)


    def add_residual(self, res_name, res):
        self.__setitem__(res_name, res)
        res.set_chain(self)
        res.set_residual_num()
    

    def keys(self):
        return self._residual_list


    def values(self):
        return [self[k] for k in self._residual_list]


    def set_pdb(self, pdb):
        self._pdb = pdb


    def provide_residual_num(self):
        self._nb_atom += 1
        return self._nb_atom

    def write(self):
        for res in self.values():
            res.write()

    def getCA(self):
        out = []
        for res in self.values():
            out += res.getCA()
        return out

    def get_chain_ID(self):
        return self._name


    def get_atom_list(self):
        ret = []
        for r in _residual_list:
            ret += r.keys()
            

class Residual(dict):
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
            value.set_residual_context(self._name, key)
        super().__setitem__(key, value)


    def add_atom(self, atomtype, atom):
        self.__setitem__(atomtype, atom)


    def keys(self):
        return self._atom_list


    def values(self):
        return [self[k] for k in self._atom_list]


    def set_chain(self, chain):
        self._chain = chain


    def set_residual_num(self):
        self._resNum = self._chain.provide_residual_num()
        

    def write(self):
        for atomtype in self.keys():
            atom = self[atomtype]
            atom.write(atomtype, self._name, self._resNum,
                       self._chain.get_chain_ID())

    def getCA(self):
        out = []
        atom = self['CA']
        out+=atom.getCoord()
        return out

        
class Atom(object):
    def __init__(self, x, y, z, identifier, symbol):
        self.x = x
        self.y = y
        self.z = z
        self._identifier = identifier
        self._residual = None
        self._symbol = symbol
        self._context = None

    
    def set_residual(self, residual):
        self._residual = residual

         
    def write(self, atomtype, resName, resSeq, chain_ID):
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



    def set_residual_context(self, res_name, atom_type):
        self._context = {"residue": res_name, "atom_type": atom_type}

    def getCoord(self):
        out = []
        out.append(self.x) # X coordiates
        out.append(self.y) # Y coordiates
        out.append(self.z) # Z coordiates
        return out


    def get_residue():
        return self._residual

class Pdb(dict):
    def __init__(self, name="pdb"):
        super().__init__()
        self._name = name
        self._chains = []
        self._nb_atom = 0


    def __setitem__(self, key, value):
        if key in self._chains:
            warn_mess = "Warning, adding " + key + " which is already in " + self._name
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


    def provide_atom_num():
        self._nb_atom += 1
        return self._nb_atom
        

class Chain(dict):
    def __init__(self, name=""):
        super().__init__()
        self._residual_list = []
        self._name = name
        self._nb_atom = 0
        self._pdb = None

    def __setitem__(self, key, value):
        if key in self._residual_list:
            warn_mess = "Warning, adding " + key + " which is already in " + self._name
            print(warn_mess, sys.stderr)
        else:
            self._residual_list.append(key)
        super().__setitem__(key, value)


    def add_residual(self, res_name, res):
        self.__setitem__(res_name, res)
        res.set_chain(self)
    

    def keys(self):
        return self._residual_list


    def values(self):
        return [self[k] for k in self._residual_list]


    def set_pdb(self, pdb):
        self._pdb = pdb


    def provide_atom_num(self):
        return self._pdb.provide_atom_num()


    def write(self):
        print("TODO")

            

class Residual(dict):
    def __init__(self, name=""):
        super().__init__()
        self._atom_list = []
        self._name = name
        self._chain = None
        self._nb_atom = 0


    def __setitem__(self, key, value):
        if key in self._atom_list:
            warn_mess = "Warning, adding " + key + " which is already in " + self._name
            print(warn_mess, sys.stderr)
        else:
            self._atom_list.append(key)
        super().__setitem__(key, value)


    def add_atom(self, atomtype, atom):
        self.__setitem__(atomtype, atom)


    def keys(self):
        return self._atom_list


    def values(self):
        return [self[k] for k in self._atom_list]


    def set_chain(self, chain):
        self._chain = chain


    def provide_atom_num(self):
        return self._chain.provide_atom_num()
        

        
class Atom(object):
    def __init__(self, x, y, z, identifier):
        self.x = x
        self.y = y
        self.z = z
        self.identifier = identifier
        self._residual = None
        self._atom_num = 0

    
    def set_residual(self, residual):
        self._residual = residual

        
    def set_atom_num(self):
        self._atom_num = self_residual.provide_atom_num()
        

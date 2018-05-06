class Chain(dict):
    def __init__(self, name=""):
        super(Chain, dict).__init__()
        self._residual_list = []
        self._name = name

    def __set_item__(self, key, value):
        if key in _residual_list:
            warn_mess = "Warning, adding " + key + " which is already in " + self._name
            print(warn_mess, sys.stderr)
        else:
            self._residual_list.append(key)
        super(Residual, dict).__set_item__(key, value)


    def add_residual(self, res_name, res):
        self.__set_item__(res_name, res)
    

    def keys(self):
        return self._residual_list


    def values(self):
        return [self[k] for k in self._residual_list]
            

class Residual(dict):
    def __init__(self, name=""):
        super(Residual, dict).__init__()
        self._atom_list = []
        self._name = name


    def __set_item__(self, key, value):
        if key in _atom_list:
            warn_mess = "Warning, adding " + key + " which is already in " + self._name
            print(warn_mess, sys.stderr)
        else:
            self._atom_list.append(key)
        super(Residual, dict).__set_item__(key, value)


    def add_atom(self, atomtype, atom):
        self.__set_item__(atomtype, atom)


    def keys(self):
        return self._atom_list


    def values(self):
        return [self[k] for k in self._atom_list]

        
class Atom(object):
    def __init__(self, x, y, z, identifier):
        self.x = x
        self.y = y
        self.z = z
        self.identifier = identifier

        

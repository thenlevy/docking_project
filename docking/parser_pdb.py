from Protein import Atom, Residual, Chain, Pdb

def parse_pdb(infile):
    """
    Reads a pdb file and return a pdb object.
    
    :param infile: the .pdb file to be read
    :type infile: str
    :return: A :py:class:`PDB` object representing the information contained
    in the .pdb file
    """

    # lecture du fichier PDB 
    f = open(infile, "r")
    lines = f.readlines()
    f.close()

    ret = Pdb()

    for line in lines:
        if line[0:4] == "ATOM":
            chain = line[21]
            if not chain in ret.keys():
                ret[chain] = Chain()
            curres = "%s"%(line[22:26]).strip()
            if not curres in ret[chain].keys():
                resname = str.strip(line[17:20])
                ret[chain].add_residual(curres, Residual(resname))
            atomtype = str.strip(line[12:16])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            identifer = line[6:11].strip()
            ret[chain][curres].add_atom(atomtype, Atom(x, y, z, identifer))

    return ret
            

from Cornell import Cornell_calc
from parser_pdb import parse_pdb
from structureToolsProjPython import preparePDB

official = preparePDB("../data/ex.pdb", "D")
test = Cornell_calc(parse_pdb("../data/ex.pdb"))

for chain in official["chains"]:
    for curres in official[chain]["reslist"]:
        for atomtype in official[chain][curres]["atomlist"]:
            atom = official[chain][curres][atomtype]
            id_at = int(atom['id']) - 1
            if atom["charge"] != test.get_q(id_at):
                print (curres, official[chain][curres]["resname"], atomtype, "charge")
                print(id_at)
            if atom["vdw"] != test.get_r(id_at):
                print (curres, official[chain][curres]["resname"], atomtype, "vdw")
                print(id_at)
            if atom["epsilon"] != test.get_epsilon(id_at):
                print (curres, official[chain][curres]["resname"], atomtype, "epsilon")
                print(id_at)

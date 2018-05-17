from parser_pdb import parse_pdb
from Protein import interface
import os


good_pdb = []
rec = parse_pdb("../data/Rec_natif.pdb")
for pdb_file in os.listdir("../data"):
    important_residues = set({'39', '42', '29', '35', '76', '38'})
    if (pdb_file.find("natif") == -1 and pdb_file != 'ex.pdb'
        and pdb_file.find(".pdb") != -1):
        print(pdb_file) # to follow the progression
        lig = parse_pdb("../data/" + pdb_file)
        if len(interface(rec.values()[0], 
                        lig.values()[0])[1].intersection(important_residues)
              ) >= 4:
            good_pdb.append(pdb_file)


good_pdb_file = open("../data/scoring_Cornell/good_pdb.txt", 'w')
for pdb in good_pdb:
    good_pdb_file.write(pdb + '\n')

        

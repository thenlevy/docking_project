from Cornell import comp_score
from RMSD import rmsdLigant
import os
from shutil import copyfile

scores = []
score_map = { }
rmsd_map = { }
good_pdb = open("../data/scoring_Cornell/good_pdb.txt", 'r')
for pdb_file in good_pdb.readlines():
    # We want to exclude the example files from the analysis
    pdb_file = pdb_file.strip()
    if (pdb_file.find("natif") == -1 and pdb_file != 'ex.pdb'
        and pdb_file.find(".pdb") != -1):
        print(pdb_file) # to follow the progression
        score = comp_score("../data/Rec_natif.pdb", "../data/" + pdb_file)
        scores.append(score)
        score_map[score] = pdb_file
        rmsd_map[score] = rmsdLigant("../data/Lig_natif.pdb", "../data/" + pdb_file)

if not os.path.exists("../data/scoring_Cornell"):
    os.mkdir("../data/scoring_Cornell")

ranking = open("../data/scoring_Cornell/ranking_filter.txt", 'w')
rank_rmsd = open("../data/scoring_Cornell/score_and_rmsd_filter.csv", 'w')

rank_rmsd.write("file;score;rmsd\n")
scores.sort()
rank = 1
for score in scores:
    ranking.write("rank {} : {} with score {}\n".format(rank, score_map[score], score))
    rank_rmsd.write("{};{};{}\n".format(score_map[score], score, rmsd_map[score]))
    rank += 1


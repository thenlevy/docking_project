from Cornell import comp_score
import os
from shutil import copyfile

scores = []
score_map = { }
for pdb_file in os.listdir("../data"):
    # We want to exclude the example files from the analysis
    if (pdb_file.find("natif") == -1 and pdb_file != 'ex.pdb'
        and pdb_file.find(".pdb") != -1):
        print(pdb_file) # to follow the progression
        score = comp_score("../data/Rec_natif.pdb", "../data/" + pdb_file)
        scores.append(score)
        score_map[score] = pdb_file

if not os.path.exists("../data/scoring_Cornell"):
    os.mkdir("../data/scoring_Cornell")

ranking = open("../data/scoring_Cornell/ranking.txt", 'w')
scores.sort()
rank = 1
for score in scores:
    ranking.write("rank {} : {} with score {}\n".format(rank, score_map[score], score))
    rank += 1

best_pdb = score_map[scores[0]]
copyfile("../data/" + best_pdb, "../data/scoring_Cornell/complexe_predit_score1.pdb")

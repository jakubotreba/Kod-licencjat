import csv
import fileinput
import IsoSpecPy
from tqdm import tqdm

def count_totalprob(mass, formula, prob):
    s=IsoSpecPy.IsoTotalProb(prob, formula)
    s.normalize()
    return (mass, formula, s)

def count_wasserstein(spect1, spect2):
    wasserstein_distance = spect1.wassersteinDistance(spect2)
    return wasserstein_distance
           
def progr(window):
    with open("../sorted_formulas.txt") as f:
        with open("results.csv", 'w') as write_file:

            writer = csv.writer(write_file)

            considered = []
            
            spectres = {}
            
            for line in tqdm(f, total = 177754527):
                mass, formula = line.split()
                mass = float(mass[:-1])
                
                if ((mass, formula)) not in considered:
                   considered.append((mass, formula))
                
                considered_filtered = [x for x in considered if x[0] > mass - window]
                
                if len(considered_filtered) > 1:
                    for x in considered_filtered:
                        isototal = count_totalprob(x[0], x[1], 0.99)
                        spectres[x[1]] = (isototal[0], isototal[2], isototal[1])
                
                for x in considered:
                    if x not in considered_filtered and x[1] in spectres:
                        del spectres[x[1]]
                                                
                considered = considered_filtered
                
                if len(considered) > 1:
                    for i in range(len(considered)):
                        for j in range(i+1, len(considered)):
                            wass = count_wasserstein(spectres[considered[i][1]][1], spectres[considered[j][1]][1])
                            mass_difference = abs(spectres[considered[i][1]][0] - spectres[considered[j][1]][0])
                            formula1, formula2 = spectres[considered[i][1]][2], spectres[considered[j][1]][2]
                            writer.writerow([mass_difference, wass, formula1, formula2])
                
progr(3.0)

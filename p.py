from random import *
from tqdm import tqdm
import linecache
import csv
import IsoSpecPy

def count_totalprob(mass, formula, prob):
    s=IsoSpecPy.IsoTotalProb(prob, formula)
    s.normalize()
    return (mass, formula, s)

def count_wasserstein(spect1, spect2):
    wasserstein_distance = spect1.wassersteinDistance(spect2)
    return wasserstein_distance
           
def progr(window):
    with open("../../sorted_formulas.txt") as f:
        with open("results205.csv", 'w') as write_file:

            writer = csv.writer(write_file)

            considered = []
            
            spectres = {}
            
            for i, line in tqdm(enumerate(f), total = 177754527):
                if random() <= 0.000005:
                    #print(i)
                    mass, formula = line.split()
                    mass = float(mass[:-1])
                    considered.append((mass, formula))
                    #considered_filtered = [x for x in considered if x[0] > mass - window]
                    counter = i
                    
                    #w tyl

                    while True:
                        line_ = linecache.getline("../../sorted_formulas.txt", counter - 1)
                        mass_, formula_ = line_.split()
                        mass_ = float(mass_[:-1])
                        if mass_ + window >= mass:
                            considered.append((mass_, formula_))
                            counter -= 1
                        else:
                            break

                    #w przod

                    while True:
                        line_ = linecache.getline("../../sorted_formulas.txt", counter + 1)
                        mass_, formula_ = line_.split()
                        mass_ = float(mass_[:-1])
                        if mass_ - window <= mass:
                            considered.append((mass_, formula_))
                            counter += 1
                        else:
                            break

                    if len(considered) > 1:
                        for x in considered:
                            isototal = count_totalprob(x[0], x[1], 0.99)
                            # mass  spectre  formula  
                            spectres[x[1]] = (isototal[0], isototal[2], isototal[1])
                        for j in range(1, len(considered)):
                            wass = count_wasserstein(spectres[considered[j][1]][1], spectres[considered[0][1]][1])
                            mass_difference = abs(spectres[considered[j][1]][0] - spectres[considered[0][1]][0])
                            formula1, formula2 = spectres[considered[j][1]][2], spectres[considered[0][1]][2]
                            writer.writerow([mass_difference, wass, formula1, formula2])
                    
                    considered = []
                    spectres = {}

                else:
                    continue            

                    """
                    for x in considered:
                        if x not in considered_filtered and x[1] in spectres:
                            del spectres[x[1]]
                                                    
                    #considered = considered_filtered

                    #the_protein = ((mass, formula))
                    
                    if len(considered) > 1:
                        for i in range(len(considered)):
                            for j in range(i+1, len(considered)):
                                wass = count_wasserstein(spectres[considered[i][1]][1], spectres[considered[j][1]][1])
                                mass_difference = abs(spectres[considered[i][1]][0] - spectres[considered[j][1]][0])
                                formula1, formula2 = spectres[considered[i][1]][2], spectres[considered[j][1]][2]
                                writer.writerow([mass_difference, wass, formula1, formula2])
                    """
progr(0.5)

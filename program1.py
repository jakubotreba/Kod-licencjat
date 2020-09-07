import fileinput
import csv
import IsoSpecPy
from tqdm import tqdm

def progr(window):
    with open('testowy.txt') as file:
        with open('test.csv', 'w', newline='') as write_file:
            writer=csv.writer(write_file)
            L=[]
            P=[]
            S=[]
            masses=[]
            formulas=[]
            limes=float('inf')
            for line in tqdm(file, total=177754527):
                average_mass=float(line.split(",")[0])
                formula=line.split(",")[1].strip()
                if average_mass<=limes:
                    L.append((average_mass, formula))
                    limes=average_mass+window
                    #print(L)
                else:
                    #rob IsoSpecPy
                    if len(L)>1:
                        #print(L)
                        combined_s=[]
                        combined_mass=[]
                        combined_formulas=[]
                        for i in range(len(L)):
                            s=IsoSpecPy.IsoTotalProb(0.99, L[i][1])
                            s.normalize()
                            mass=L[i][0]
                            S.append(s)
                            masses.append(mass)
                            formulas.append(L[i][1])
                        for i in range(len(L)):
                            for j in range(i+1, len(L)):
                                combined_s.append((S[i], S[j]))
                                combined_mass.append((masses[i], masses[j]))
                                combined_formulas.append((formulas[i], formulas[j]))
                        S=[]
                        masses=[]
                        #print(len(combined_s))
                        for i in range(len(combined_s)):
                            wasserstein=combined_s[i][0].wassersteinDistance(combined_s[i][1])
                            mass_difference=abs(combined_mass[i][0]-combined_mass[i][1])
                            used_formulas=combined_formulas[i]
                            #print(mass_difference, wasserstein, used_formulas)
                            writer.writerow([(mass_difference, wasserstein, used_formulas)])
                            #print("napisalem")
                        formulas=[]
                        combined_s=[]
                        combined_mass=[]
                        combined_formulas=[]
                    L.append((average_mass, formula))
                    for i in range(len(L)):
                        if L[i][0]==average_mass:
                            P.append(L[i])
                            #print(P)
                    L=P
                    P=[]
                    limes=average_mass+window
                    #print(L)
progr(10.0)

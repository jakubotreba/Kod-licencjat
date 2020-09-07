from tqdm import tqdm
from random import *
import IsoSpecPy
import IsoSpecPy.Distributions
import numpy as np

def grouper():
	
	with open("../../sorted_formulas.txt") as f:
		prev = None
		group = []
		for item in tqdm(f, total = 150419758):
			if not prev or float(item.split()[0].strip(",")) <= float(group[0][0]) + 0.5:
				item1 = item.split(",")[0].strip(",")
				item2 = item.split(",")[1].strip()
				group.append([item1, item2])
			else:
				yield group
				group = [[item.split(",")[0].strip(","), item.split(",")[1].strip()]]
			prev = item
		if group:
			yield group

def stochastic(clusters, mol_number, cluster_n, iteration):
	with open(str(iteration) + "stochastic" + str(mol_number) + ".txt", "w+") as g:
		for cluster in [clusters[cluster_n]]:
			the_protein = cluster[int(len(cluster) / 2)]
			the_protein_name = the_protein[1]
			print(len(cluster), "stochastic")
			A = IsoSpecPy.IsoTotalProb(0.99, the_protein_name)
			A.normalize()
			Astoch = IsoSpecPy.IsoStochasticGenerator(mol_number, the_protein_name)
			M, P = zip(*Astoch)
			W = IsoSpecPy.IsoDistribution(masses = M, probs = P)
			W.normalize()
			g.write("<" + str(A.wassersteinDistance(W)) + ", " + the_protein_name + ", " + the_protein_name + "\n")
			for protein in cluster:
				if protein == the_protein:
					continue
				protein_name = protein[1]
				a = IsoSpecPy.IsoTotalProb(0.99, protein_name)
				a.normalize()
				g.write(str(a.wassersteinDistance(W)) + ", " + protein_name + "; " + the_protein_name + "\n")
			break
	g.close()

def gaussian(clusters, stdev, cluster_n, iteration):
	
	with open(str(iteration) + "gaussians" + str(stdev) + ".txt", "w+") as g:
		G = IsoSpecPy.Distributions.Gaussian(stdev, 0.001, 0.99)
		for cluster in [clusters[cluster_n]]:
			print(len(cluster), "gaussian")
			the_protein = cluster[int(len(cluster) / 2)]
			the_protein_name = the_protein[1]
			the_protein_mass = the_protein[0]
			A = IsoSpecPy.IsoTotalProb(0.99, the_protein_name)
			AB = A.binned(0.01)
			AB_unsettled = AB * G
			g.write("<" + str(A.wassersteinDistance(AB_unsettled)) + ", " + the_protein_name + ", " + the_protein_name + "\n")
			for protein in cluster:
				if protein == the_protein:
					continue
				protein_name = protein[1]
				protein_mass = protein[0]
				a = IsoSpecPy.IsoTotalProb(0.99, protein_name)
				g.write(str(a.wassersteinDistance(AB_unsettled)) + ", " + protein_name + "; " + the_protein_name + "\n")
			break        
	g.close()

def electronic(clusters, signal, noise, cluster_n, iteration):
	
	with open(str(iteration) + "electronic" + str(signal) + ".txt", "w+") as g:
		for cluster in [clusters[cluster_n]]:
			the_protein = cluster[int(len(cluster) / 2)]
			the_protein_name = the_protein[1]
			print(len(cluster), "electronic")
			A = IsoSpecPy.IsoTotalProb(0.99, the_protein_name)
			AB = A.binned(0.1)
			M = []
			for mass, prob in AB:
				M.append((mass, prob))
			M.sort()
			min_mass = M[0][0]
			max_mass = M[-1][0]
			masses = list(np.arange(min_mass - 0.1, max_mass + 0.1, 0.1))
			probs = [np.random.poisson(100.0) for _ in masses]
			S = IsoSpecPy.IsoDistribution(masses = masses, probs = probs)
			S.normalize()
			AB_unsettled = IsoSpecPy.IsoDistribution.LinearCombination([A, S], [signal, noise])
			AB_uns = AB_unsettled.binned(0.1)
			AB_uns.normalize()
			A.normalize()
			g.write("<" + str(A.wassersteinDistance(AB_uns)) + ", " + the_protein_name + ", " + the_protein_name + "\n")
			for protein in cluster:
				if protein == the_protein:
					continue
				protein_name = protein[1]
				a = IsoSpecPy.IsoTotalProb(0.99, protein_name)
				a.normalize()
				g.write(str(a.wassersteinDistance(AB_uns)) + ", " + protein_name + "; " + the_protein_name + "\n")
			break        
	g.close()

def sorting_clusters():
	
	print("rozpoczalem sortowanie klastrow")
	L = []
	D = dict(enumerate(grouper(), 1))
	for k in sorted(D, key = lambda k: len(D[k]), reverse = True):
		L.append(D[k])

	print("skonczylem sortowanie klastrow")
	return L

klastry = sorting_clusters()
print(len(klastry[0]), klastry[0][0])
print(len(klastry[4000]), klastry[4000][0])
print(len(klastry[81000]), klastry[81000][0])
print(len(klastry[90000]), klastry[90000][0])
print(len(klastry[189000]), klastry[189000][0])
print(len(klastry[190000]), klastry[190000][0])

#wektor = np.linspace(0, 6, num = 10)
#os_x = [10**w for w in wektor]

"""
#1
for i in os_x:
	stochastic(klastry, int(i), 0, 1)
#2
for i in os_x:
	stochastic(klastry, int(i), 4000, 2)
#3
for i in os_x:
	stochastic(klastry, int(i), 81000, 3)
#4
for i in os_x:
	stochastic(klastry, int(i), 90000, 4)
#5
for i in os_x:
	stochastic(klastry, int(i), 189000, 5)
#6
for i in os_x:
	stochastic(klastry, int(i), 190000, 6)
"""
"""
1. [,], 1750 0 
2.['9447.200036444021', 'C360H691N99O182S4'] 1510 4000
3.['46884.5343400998', 'C1741H3546N530O920S8'] 860 81000
4.['8625.868884764279', 'C352H676N78O156S4'] 755 90000
5.['100417.19805292874', 'C3820H7596N1044O1964S26'] 135 189000
6. [,] ? 190000
"""


#1
electronic(klastry, 1.0, 0.0, 0, 1)
electronic(klastry, 0.95, 0.05, 0, 1)
electronic(klastry, 0.9, 0.1, 0, 1)
electronic(klastry, 0.85, 0.15, 0, 1)
electronic(klastry, 0.80, 0.20, 0, 1)
electronic(klastry, 0.75, 0.25, 0, 1)
electronic(klastry, 0.70, 0.30, 0, 1)
electronic(klastry, 0.65, 0.35, 0, 1)
electronic(klastry, 0.60, 0.40, 0, 1)
electronic(klastry, 0.55, 0.45, 0, 1)
electronic(klastry, 0.50, 0.50, 0, 1)
#2
electronic(klastry, 1.0, 0.0, 4000, 2)
electronic(klastry, 0.95, 0.05, 4000, 2)
electronic(klastry, 0.9, 0.1, 4000, 2)
electronic(klastry, 0.85, 0.15, 4000, 2)
electronic(klastry, 0.80, 0.20, 4000, 2)
electronic(klastry, 0.75, 0.25, 4000, 2)
electronic(klastry, 0.70, 0.30, 4000, 2)
electronic(klastry, 0.65, 0.35, 4000, 2)
electronic(klastry, 0.60, 0.40, 4000, 2)
electronic(klastry, 0.55, 0.45, 4000, 2)
electronic(klastry, 0.50, 0.50, 4000, 2)
#3
electronic(klastry, 1.0, 0.0, 81000, 3)
electronic(klastry, 0.95, 0.05, 81000, 3)
electronic(klastry, 0.9, 0.1, 81000, 3)
electronic(klastry, 0.85, 0.15, 81000, 3)
electronic(klastry, 0.80, 0.20, 81000, 3)
electronic(klastry, 0.75, 0.25, 81000, 3)
electronic(klastry, 0.70, 0.30, 81000, 3)
electronic(klastry, 0.65, 0.35, 81000, 3)
electronic(klastry, 0.60, 0.40, 81000, 3)
electronic(klastry, 0.55, 0.45, 81000, 3)
electronic(klastry, 0.50, 0.50, 81000, 3)
#4
electronic(klastry, 1.0, 0.0, 90000, 4)
electronic(klastry, 0.95, 0.05, 90000, 4)
electronic(klastry, 0.9, 0.1, 90000, 4)
electronic(klastry, 0.85, 0.15, 90000, 4)
electronic(klastry, 0.80, 0.20, 90000, 4)
electronic(klastry, 0.75, 0.25, 90000, 4)
electronic(klastry, 0.70, 0.30, 90000, 4)
electronic(klastry, 0.65, 0.35, 90000, 4)
electronic(klastry, 0.60, 0.40, 90000, 4)
electronic(klastry, 0.55, 0.45, 90000, 4)
electronic(klastry, 0.50, 0.50, 90000, 4)
#5
electronic(klastry, 1.0, 0.0, 189000, 5)
electronic(klastry, 0.95, 0.05, 189000, 5)
electronic(klastry, 0.9, 0.1, 189000, 5)
electronic(klastry, 0.85, 0.15, 189000, 5)
electronic(klastry, 0.80, 0.20, 189000, 5)
electronic(klastry, 0.75, 0.25, 189000, 5)
electronic(klastry, 0.70, 0.30, 189000, 5)
electronic(klastry, 0.65, 0.35, 189000, 5)
electronic(klastry, 0.60, 0.40, 189000, 5)
electronic(klastry, 0.55, 0.45, 189000, 5)
electronic(klastry, 0.50, 0.50, 189000, 5)
#6
electronic(klastry, 1.0, 0.0, 190000, 6)
electronic(klastry, 0.95, 0.05, 190000, 6)
electronic(klastry, 0.9, 0.1, 190000, 6)
electronic(klastry, 0.85, 0.15, 190000, 6)
electronic(klastry, 0.80, 0.20, 190000, 6)
electronic(klastry, 0.75, 0.25, 190000, 6)
electronic(klastry, 0.70, 0.30, 190000, 6)
electronic(klastry, 0.65, 0.35, 190000, 6)
electronic(klastry, 0.60, 0.40, 190000, 6)
electronic(klastry, 0.55, 0.45, 190000, 6)
electronic(klastry, 0.50, 0.50, 190000, 6)


"""
#1
gaussian(klastry, 1e-06, 0, 1)
gaussian(klastry, 1.778279410038923e-06, 0, 1)
gaussian(klastry, 3.162277660168379e-06, 0, 1)
gaussian(klastry, 5.623413251903491e-06, 0, 1)
gaussian(klastry, 1e-05, 0, 1)
gaussian(klastry, 1.778279410038923e-05, 0, 1)
gaussian(klastry, 3.1622776601683795e-05, 0, 1)
gaussian(klastry, 5.623413251903491e-05, 0, 1)
gaussian(klastry, 0.0001, 0, 1)
gaussian(klastry, 0.00017782794100389227, 0, 1)
gaussian(klastry, 0.00031622776601683794, 0, 1)
gaussian(klastry, 0.0005623413251903491, 0, 1)
gaussian(klastry, 0.001, 0, 1)
gaussian(klastry, 0.0017782794100389228, 0, 1)
gaussian(klastry, 0.0031622776601683794, 0, 1)
gaussian(klastry, 0.005623413251903491, 0, 1)
gaussian(klastry, 0.01, 0, 1)
gaussian(klastry, 0.01778279410038923, 0, 1)
gaussian(klastry, 0.03162277660168379, 0, 1)
gaussian(klastry, 0.05623413251903491, 0, 1)
gaussian(klastry, 0.1, 0, 1)
gaussian(klastry, 0.1778279410038923, 0, 1)
gaussian(klastry, 0.31622776601683794, 0, 1)
gaussian(klastry, 0.5623413251903491, 0, 1)
gaussian(klastry, 1.0, 0, 1)
gaussian(klastry, 1.7782794100389228, 0, 1)
gaussian(klastry, 3.1622776601683795, 0, 1)
gaussian(klastry, 5.623413251903491, 0, 1)
gaussian(klastry, 10.0, 0, 1)
#2
gaussian(klastry, 1e-06, 4000, 2)
gaussian(klastry, 1.778279410038923e-06, 4000, 2)
gaussian(klastry, 3.162277660168379e-06, 4000, 2)
gaussian(klastry, 5.623413251903491e-06, 4000, 2)
gaussian(klastry, 1e-05, 4000, 2)
gaussian(klastry, 1.778279410038923e-05, 4000, 2)
gaussian(klastry, 3.1622776601683795e-05, 4000, 2)
gaussian(klastry, 5.623413251903491e-05, 4000, 2)
gaussian(klastry, 0.0001, 4000, 2)
gaussian(klastry, 0.00017782794100389227, 4000, 2)
gaussian(klastry, 0.00031622776601683794, 4000, 2)
gaussian(klastry, 0.0005623413251903491, 4000, 2)
gaussian(klastry, 0.001, 4000, 2)
gaussian(klastry, 0.0017782794100389228, 4000, 2)
gaussian(klastry, 0.0031622776601683794, 4000, 2)
gaussian(klastry, 0.005623413251903491, 4000, 2)
gaussian(klastry, 0.01, 4000, 2)
gaussian(klastry, 0.01778279410038923, 4000, 2)
gaussian(klastry, 0.03162277660168379, 4000, 2)
gaussian(klastry, 0.05623413251903491, 4000, 2)
gaussian(klastry, 0.1, 4000, 2)
gaussian(klastry, 0.1778279410038923, 4000, 2)
gaussian(klastry, 0.31622776601683794, 4000, 2)
gaussian(klastry, 0.5623413251903491, 4000, 2)
gaussian(klastry, 1.0, 4000, 2)
gaussian(klastry, 1.7782794100389228, 4000, 2)
gaussian(klastry, 3.1622776601683795, 4000, 2)
gaussian(klastry, 5.623413251903491, 4000, 2)
gaussian(klastry, 10.0, 4000, 2)
#3
gaussian(klastry, 1e-06, 81000, 3)
gaussian(klastry, 1.778279410038923e-06, 81000, 3)
gaussian(klastry, 3.162277660168379e-06, 81000, 3)
gaussian(klastry, 5.623413251903491e-06, 81000, 3)
gaussian(klastry, 1e-05, 81000, 3)
gaussian(klastry, 1.778279410038923e-05, 81000, 3)
gaussian(klastry, 3.1622776601683795e-05, 81000, 3)
gaussian(klastry, 5.623413251903491e-05, 81000, 3)
gaussian(klastry, 0.0001, 81000, 3)
gaussian(klastry, 0.00017782794100389227, 81000, 3)
gaussian(klastry, 0.00031622776601683794, 81000, 3)
gaussian(klastry, 0.0005623413251903491, 81000, 3)
gaussian(klastry, 0.001, 81000, 3)
gaussian(klastry, 0.0017782794100389228, 81000, 3)
gaussian(klastry, 0.0031622776601683794, 81000, 3)
gaussian(klastry, 0.005623413251903491, 81000, 3)
gaussian(klastry, 0.01, 81000, 3)
gaussian(klastry, 0.01778279410038923, 81000, 3)
gaussian(klastry, 0.03162277660168379, 81000, 3)
gaussian(klastry, 0.05623413251903491, 81000, 3)
gaussian(klastry, 0.1, 81000, 3)
gaussian(klastry, 0.1778279410038923, 81000, 3)
gaussian(klastry, 0.31622776601683794, 81000, 3)
gaussian(klastry, 0.5623413251903491, 81000, 3)
gaussian(klastry, 1.0, 81000, 3)
gaussian(klastry, 1.7782794100389228, 81000, 3)
gaussian(klastry, 3.1622776601683795, 81000, 3)
gaussian(klastry, 5.623413251903491, 81000, 3)
gaussian(klastry, 10.0, 81000, 3)
#4
gaussian(klastry, 1e-06, 90000, 4)
gaussian(klastry, 1.778279410038923e-06, 90000, 4)
gaussian(klastry, 3.162277660168379e-06, 90000, 4)
gaussian(klastry, 5.623413251903491e-06, 90000, 4)
gaussian(klastry, 1e-05, 90000, 4)
gaussian(klastry, 1.778279410038923e-05, 90000, 4)
gaussian(klastry, 3.1622776601683795e-05, 90000, 4)
gaussian(klastry, 5.623413251903491e-05, 90000, 4)
gaussian(klastry, 0.0001, 90000, 4)
gaussian(klastry, 0.00017782794100389227, 90000, 4)
gaussian(klastry, 0.00031622776601683794, 90000, 4)
gaussian(klastry, 0.0005623413251903491, 90000, 4)
gaussian(klastry, 0.001, 90000, 4)
gaussian(klastry, 0.0017782794100389228, 90000, 4)
gaussian(klastry, 0.0031622776601683794, 90000, 4)
gaussian(klastry, 0.005623413251903491, 90000, 4)
gaussian(klastry, 0.01, 90000, 4)
gaussian(klastry, 0.01778279410038923, 90000, 4)
gaussian(klastry, 0.03162277660168379, 90000, 4)
gaussian(klastry, 0.05623413251903491, 90000, 4)
gaussian(klastry, 0.1, 90000, 4)
gaussian(klastry, 0.1778279410038923, 90000, 4)
gaussian(klastry, 0.31622776601683794, 90000, 4)
gaussian(klastry, 0.5623413251903491, 90000, 4)
gaussian(klastry, 1.0, 90000, 4)
gaussian(klastry, 1.7782794100389228, 90000, 4)
gaussian(klastry, 3.1622776601683795, 90000, 4)
gaussian(klastry, 5.623413251903491, 90000, 4)
gaussian(klastry, 10.0, 90000, 4)
#5
gaussian(klastry, 1e-06, 189000, 5)
gaussian(klastry, 1.778279410038923e-06, 189000, 5)
gaussian(klastry, 3.162277660168379e-06, 189000, 5)
gaussian(klastry, 5.623413251903491e-06, 189000, 5)
gaussian(klastry, 1e-05, 189000, 5)
gaussian(klastry, 1.778279410038923e-05, 189000, 5)
gaussian(klastry, 3.1622776601683795e-05, 189000, 5)
gaussian(klastry, 5.623413251903491e-05, 189000, 5)
gaussian(klastry, 0.0001, 189000, 5)
gaussian(klastry, 0.00017782794100389227, 189000, 5)
gaussian(klastry, 0.00031622776601683794, 189000, 5)
gaussian(klastry, 0.0005623413251903491, 189000, 5)
gaussian(klastry, 0.001, 189000, 5)
gaussian(klastry, 0.0017782794100389228, 189000, 5)
gaussian(klastry, 0.0031622776601683794, 189000, 5)
gaussian(klastry, 0.005623413251903491, 189000, 5)
gaussian(klastry, 0.01, 189000, 5)
gaussian(klastry, 0.01778279410038923, 189000, 5)
gaussian(klastry, 0.03162277660168379, 189000, 5)
gaussian(klastry, 0.05623413251903491, 189000, 5)
gaussian(klastry, 0.1, 189000, 5)
gaussian(klastry, 0.1778279410038923, 189000, 5)
gaussian(klastry, 0.31622776601683794, 189000, 5)
gaussian(klastry, 0.5623413251903491, 189000, 5)
gaussian(klastry, 1.0, 189000, 5)
gaussian(klastry, 1.7782794100389228, 189000, 5)
gaussian(klastry, 3.1622776601683795, 189000, 5)
gaussian(klastry, 5.623413251903491, 189000, 5)
gaussian(klastry, 10.0, 189000, 5)
#6
gaussian(klastry, 1e-06, 190000, 6)
gaussian(klastry, 1.778279410038923e-06, 190000, 6)
gaussian(klastry, 3.162277660168379e-06, 190000, 6)
gaussian(klastry, 5.623413251903491e-06, 190000, 6)
gaussian(klastry, 1e-05, 190000, 6)
gaussian(klastry, 1.778279410038923e-05, 190000, 6)
gaussian(klastry, 3.1622776601683795e-05, 190000, 6)
gaussian(klastry, 5.623413251903491e-05, 190000, 6)
gaussian(klastry, 0.0001, 190000, 6)
gaussian(klastry, 0.00017782794100389227, 190000, 6)
gaussian(klastry, 0.00031622776601683794, 190000, 6)
gaussian(klastry, 0.0005623413251903491, 190000, 6)
gaussian(klastry, 0.001, 190000, 6)
gaussian(klastry, 0.0017782794100389228, 190000, 6)
gaussian(klastry, 0.0031622776601683794, 190000, 6)
gaussian(klastry, 0.005623413251903491, 190000, 6)
gaussian(klastry, 0.01, 190000, 6)
gaussian(klastry, 0.01778279410038923, 190000, 6)
gaussian(klastry, 0.03162277660168379, 190000, 6)
gaussian(klastry, 0.05623413251903491, 190000, 6)
gaussian(klastry, 0.1, 190000, 6)
gaussian(klastry, 0.1778279410038923, 190000, 6)
gaussian(klastry, 0.31622776601683794, 190000, 6)
gaussian(klastry, 0.5623413251903491, 190000, 6)
gaussian(klastry, 1.0, 190000, 6)
gaussian(klastry, 1.7782794100389228, 190000, 6)
gaussian(klastry, 3.1622776601683795, 190000, 6)
gaussian(klastry, 5.623413251903491, 190000, 6)
gaussian(klastry, 10.0, 190000, 6)
"""

import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd 

def get_gaussian_data(files):
	
	stdevs = []
	ys = []
	
	for file in files:
		stdev0 = file[10:]
		stdev = stdev0.split("txt")[0].strip(".")
		stdevs.append(stdev)

	nr_klastra = files[0][0]

	for i in range(len(stdevs)):
		if "-" in stdevs[i]:
			if "." not in stdevs[i]:
				stdevs[i] = stdevs[i][0] + ".0" + stdevs[i][1:]
			number = stdevs[i]
			lead, power = number.split("e-")
			a, b = lead.split(".")
			number = "0."+"0"*(int(power) - 1) + a + b
			stdevs[i] = float(number)
		else:
			stdevs[i] = float(stdevs[i])

	for i in range(len(files)):
		with open(files[i]) as f:
			y = []
			for line in f:
				if "<" in line:
					y.append(float(line.split(",")[0].strip("<")))
				else:
					y.append(float(line.split(",")[0]))
		
		f.close()
		ys.append(y)

	df = pd.DataFrame(ys, index = stdevs)
	df_main = df[0]
	df = df.transpose()
	df = df.sort_values(by = stdevs[1])

	fig, ax = plt.subplots()
	
	plt.plot(df.transpose()[int(0.1*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.2*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.3*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.4*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.5*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.6*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.7*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.8*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.9*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(1.0*len(df)) - 1], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df_main.transpose(), color = 'red', linewidth=1.75)
	plt.xscale('log')
	plt.title(f'Klaster nr {nr_klastra} - szum Gaussowski')
	plt.xlabel('Odchylenie standardowe (Da)')
	plt.ylabel('Odległość Wassersteina')
	plt.savefig(f'gauss{nr_klastra}.png')

def get_electronic_data(files):
	
	noises = []
	ys = []

	for file in files:
		nois0 = file[11:]
		nois = nois0.split("txt")[0].strip(".")
		noises.append(1.0 - float(nois))

	nr_klastra = files[0][0]

	for i in range(len(files)):
		with open(files[i]) as f:
			y = []
			for line in f:
				if "<" in line:
					y.append(float(line.split(",")[0].strip("<")))
				else:
					y.append(float(line.split(",")[0]))
		
		f.close()
		ys.append(y)

	df = pd.DataFrame(ys, index = noises)
	df_main = df[0]
	df = df.transpose()
	df = df.sort_values(by = noises[1])

	fig, ax = plt.subplots()

	plt.plot(df.transpose()[int(0.1*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.2*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.3*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.4*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.5*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.6*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.7*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.8*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.9*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(1.0*len(df)) - 1], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df_main.transpose(), color = 'red', linewidth=1.75)
	plt.title(f'Klaster nr {nr_klastra} - szum elektroniczny')
	plt.xlabel('Poziom szumu')
	plt.ylabel('Odległość Wassersteina')
	plt.savefig(f'electronic{nr_klastra}.png')


def get_stochastic_data(files):
	
	molecules = []
	ys = []

	for file in files:
		molecules0 = file[11:]
		molecule = molecules0.split("txt")[0].strip(".")
		molecules.append(float(molecule))

	nr_klastra = files[0][0]

	for i in range(len(files)):
		with open(files[i]) as f:
			y = []
			for line in f:
				if "<" in line:
					y.append(float(line.split(",")[0].strip("<")))
				else:
					y.append(float(line.split(",")[0]))
		
		f.close()
		ys.append(y)

	df = pd.DataFrame(ys, index = molecules)
	df_main = df[0]
	df = df.transpose()
	df = df.sort_values(by = molecules[1])

	fig, ax = plt.subplots()

	plt.plot(df.transpose()[int(0.1*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.2*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.3*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.4*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.5*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.6*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.7*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.8*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(0.9*len(df))], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df.transpose()[int(1.0*len(df)) - 1], alpha = 0.5, linewidth=0.5, color = 'blue')
	plt.plot(df_main.transpose(), color = 'red', linewidth=1.75)
	plt.title(f'Klaster nr {nr_klastra} - szum stochastyczny')
	plt.xscale('log')
	plt.xlabel('Liczba molekuł (log)')
	plt.ylabel('Odległość Wassersteina')
	plt.savefig(f'stochastic{nr_klastra}.png')

#electronic
for i in range(1, 7):
	get_electronic_data([f"{i}electronic1.0.txt" ,f"{i}electronic0.95.txt" ,f"{i}electronic0.9.txt", f"{i}electronic0.85.txt", f"{i}electronic0.8.txt", f"{i}electronic0.75.txt", f"{i}electronic0.7.txt", f"{i}electronic0.65.txt", f"{i}electronic0.6.txt", f"{i}electronic0.55.txt", f"{i}electronic0.5.txt"])

#gaussian
#for i in range(1, 7):
#	get_gaussian_data([f"{i}gaussians1e-06.txt", f"{i}gaussians1.778279410038923e-06.txt", f"{i}gaussians3.162277660168379e-06.txt", f"{i}gaussians5.623413251903491e-06.txt", f"{i}gaussians1e-05.txt", f"{i}gaussians1.778279410038923e-05.txt", f"{i}gaussians3.1622776601683795e-05.txt", f"{i}gaussians5.623413251903491e-05.txt", f"{i}gaussians0.0001.txt", f"{i}gaussians0.00017782794100389227.txt", f"{i}gaussians0.00031622776601683794.txt", f"{i}gaussians0.0005623413251903491.txt", f"{i}gaussians0.001.txt", f"{i}gaussians0.0017782794100389228.txt", f"{i}gaussians0.0031622776601683794.txt", f"{i}gaussians0.005623413251903491.txt", f"{i}gaussians0.01.txt", f"{i}gaussians0.01778279410038923.txt", f"{i}gaussians0.03162277660168379.txt", f"{i}gaussians0.05623413251903491.txt", f"{i}gaussians0.1.txt", f"{i}gaussians0.1778279410038923.txt", f"{i}gaussians0.31622776601683794.txt", f"{i}gaussians0.5623413251903491.txt", f"{i}gaussians1.0.txt", f"{i}gaussians1.7782794100389228.txt", f"{i}gaussians3.1622776601683795.txt", f"{i}gaussians5.623413251903491.txt", f"{i}gaussians10.0.txt"])	

#stochastic
#for i in range(1, 7):
#	get_stochastic_data([f"{i}stochastic1.txt", f"{i}stochastic4.txt", f"{i}stochastic21.txt", f"{i}stochastic100.txt", f"{i}stochastic464.txt", f"{i}stochastic2154.txt", f"{i}stochastic10000.txt", f"{i}stochastic46415.txt", f"{i}stochastic215443.txt", f"{i}stochastic1000000.txt"])






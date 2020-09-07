from tqdm import tqdm
from random import *
import IsoSpecPy
import IsoSpecPy.Distributions
"""
def binning():
	with open("clusters.txt") as f:
		with open("gaussians.txt", "w+") as g:
			G = IsoSpecPy.Distributions.Gaussian(0.01, 0.001, 0.99)
			for line in tqdm(f, total = 7190607478):
				#print(line.strip())
				cluster = 
				#print(list(cluster.split(",")))
				the_protein = cluster.split()[int(len(cluster.split()) / 2)]
				#print(the_protein)
				the_protein_name = the_protein[1]
				the_protein_mass = the_protein[0]
				A = IsoSpecPy.IsoTotalProb(0.99, the_protein_name)
				AB = A.binned(0.01)
				AB_unsettled = AB * G
				for protein in cluster:
					if protein == the_protein:
						continue
					protein_name = protein[1]
					protein_mass = protein[0]
					a = IsoSpecPy.IsoTotalProb(0.99, protein_name)
					g.write(str(a.wassersteinDistance(AB_unsettled)) + ", " + protein_name + "; " + the_protein_name + "\n")
		g.close()
	f.close()

binning()
"""

L = []

with open("clusters.txt") as f:
	if 'C644H1339N195O311S7' in f:
		print("Ubikwityna!")

"""
formula=[0,0,0,0,0]
record_seq = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGIIEPSLRQLAQKYNCDKMICRKCYARLHPRAVNCRKKKCGHTNNLRPKKKVK'
for j in range(len(record_seq)):
	if record_seq[j]=='A':
		formula[0]+=3
		formula[1]+=7
		formula[2]+=1
		formula[3]+=2
	elif record_seq[j]=='C':
		formula[0]+=3
		formula[1]+=7
		formula[2]+=1
		formula[3]+=2
		formula[4]+=1
	elif record_seq[j]=='D':
		formula[0]+=4
		formula[1]+=7
		formula[2]+=1
		formula[3]+=4
	elif record_seq[j]=='E':
		formula[0]+=5
		formula[1]+=9
		formula[2]+=1
		formula[3]+=4
	elif record_seq[j]=='F':
		formula[0]+=9
		formula[1]+=11
		formula[2]+=1
		formula[3]+=2
	elif record_seq[j]=='G':
		formula[0]+=2
		formula[1]+=5
		formula[2]+=1
		formula[3]+=2
	elif record_seq[j]=='H':
		formula[0]+=6
		formula[1]+=9
		formula[2]+=3
		formula[3]+=2
	elif record_seq[j]=='I':
		formula[0]+=6
		formula[1]+=13
		formula[2]+=1
		formula[3]+=2
	elif record_seq[j]=='K':
		formula[0]+=6
		formula[1]+=14
		formula[2]+=2
		formula[3]+=2
	elif record_seq[j]=='L':
		formula[0]+=6
		formula[1]+=13
		formula[2]+=1
		formula[3]+=2
	elif record_seq[j]=='M':
		formula[0]+=5
		formula[1]+=11
		formula[2]+=1
		formula[3]+=2
		formula[4]+=1
	elif record_seq[j]=='N':
		formula[0]+=4
		formula[1]+=8
		formula[2]+=2
		formula[3]+=3
	elif record_seq[j]=='P':
		formula[0]+=5
		formula[1]+=9
		formula[2]+=1
		formula[3]+=2
	elif record_seq[j]=='Q':
		formula[0]+=5
		formula[1]+=10
		formula[2]+=2
		formula[3]+=3
	elif record_seq[j]=='R':
		formula[0]+=6
		formula[1]+=14
		formula[2]+=4
		formula[3]+=2
	elif record_seq[j]=='S':
		formula[0]+=3
		formula[1]+=7
		formula[2]+=1
		formula[3]+=3
	elif record_seq[j]=='T':
		formula[0]+=4
		formula[1]+=9
		formula[2]+=1
		formula[3]+=3
	elif record_seq[j]=='V':
		formula[0]+=5
		formula[1]+=11
		formula[2]+=1
		formula[3]+=2
	elif record_seq[j]=='W':
		formula[0]+=11
		formula[1]+=12
		formula[2]+=2
		formula[3]+=2
	elif record_seq[j]=='Y':
		formula[0]+=9
		formula[1]+=11
		formula[2]+=1
		formula[3]+=3

formula[0]='C'+str(formula[0])
formula[1]='H'+str(formula[1])
formula[2]='N'+str(formula[2])
formula[3]='O'+str(formula[3])
formula[4]='S'+str(formula[4])
formula=''.join(formula)
print(formula)
"""



















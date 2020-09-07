import IsoSpecPy
from tqdm import tqdm


L = []
seen = set()
with open("results_fasta_to_formula.txt") as f:
    for line in tqdm(f, total = 355509054):
        if ">" in line:
            continue
        if line in seen: 
            continue
        seen.add(line)
        formula = line.strip()
	
        mass = IsoSpecPy.Iso(formula).getTheoreticalAverageMass()
        L.append((mass, formula))


L.sort()
with open('sorted_formulas.txt', 'w') as f2:
    for item in L:
       f2.write(str(item[0])+', '+str(item[1])+'\n')
f2.close()

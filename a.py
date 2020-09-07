from tqdm import tqdm

def opening_clusters():
	lista = []
	with open("clusters.txt") as c:
		for item in tqdm(c):
			e = item.split()
			print(e[0])
			break

opening_clusters()
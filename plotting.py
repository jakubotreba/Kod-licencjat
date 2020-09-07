import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd 
import numpy as np

def get_data(files):
	
	stdevs = []
	ys = []
	
	for file in files:
		stdev0 = file[9:]
		stdev = stdev0.split("txt")[0].strip(".")
		stdevs.append(stdev)

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
	plt.title('Klaster nr 5 - szum Gaussowski')
	plt.xlabel('Odchylenie standardowe')
	plt.ylabel('Odległość Wassersteina')
	plt.savefig('gauss5.png')

def get_electronic_data(files):
	noises = []
	ys = []

	for file in files:
		nois0 = file[10:]
		nois = nois0.split("txt")[0].strip(".")
		noises.append(1.0 - float(nois))

	print(noises)

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
	plt.title('Klaster nr 5 - szum elektroniczny')
	plt.xlabel('Procent szumu')
	plt.ylabel('Odległość Wassersteina')
	plt.savefig('electronic5.png')


get_electronic_data(["electronic1.0.txt" ,"electronic0.95.txt" ,"electronic0.9.txt", "electronic0.85.txt", "electronic0.8.txt", "electronic0.75.txt", "electronic0.7.txt", "electronic0.65.txt", "electronic0.6.txt", "electronic0.55.txt", "electronic0.5.txt"])




"""	
	for i in range(len(ys[0])):
		for j in range(len(stdevs)):
			if i == 0:
				main_x.append(stdevs[j])
				main_y.append(ys[j][i])
			else:	
				os_x.append(stdevs[j])
				os_y.append(ys[j][i])
	
	df = pd.DataFrame(os_y)
	#df['Decile_rank'] = pd.qcut(df, 10, labels = False)

	plt.plot(os_x, os_y, alpha = 0.5)
	plt.plot(main_x, main_y, color = "red")
	plt.xscale('log')
	plt.savefig("plotGauss.png")		
"""
#get_data(["gaussians1e-06.txt", "gaussians1.778279410038923e-06.txt", "gaussians3.162277660168379e-06.txt", "gaussians5.623413251903491e-06.txt", "gaussians1e-05.txt", "gaussians1.778279410038923e-05.txt", "gaussians3.1622776601683795e-05.txt", "gaussians5.623413251903491e-05.txt", "gaussians0.0001.txt", "gaussians0.00017782794100389227.txt", "gaussians0.00031622776601683794.txt", "gaussians0.0005623413251903491.txt", "gaussians0.001.txt", "gaussians0.0017782794100389228.txt", "gaussians0.0031622776601683794.txt", "gaussians0.005623413251903491.txt", "gaussians0.01.txt", "gaussians0.01778279410038923.txt", "gaussians0.03162277660168379.txt", "gaussians0.05623413251903491.txt", "gaussians0.1.txt", "gaussians0.1778279410038923.txt", "gaussians0.31622776601683794.txt", "gaussians0.5623413251903491.txt", "gaussians1.0.txt", "gaussians1.7782794100389228.txt", "gaussians3.1622776601683795.txt", "gaussians5.623413251903491.txt", "gaussians10.0.txt"])	

"""
#stdevs = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0]

fig, ax = plt.subplots()
plot = ax.scatter(stdevs, wassersteins)
plt.savefig("plotclustertest.png")
"""
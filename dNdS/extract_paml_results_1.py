'''
Purpose: extract dN/dS ratios from codeml results, writes file with dN/dS differences
Usage: extract_paml_results.py codeml_results.txt
	Where codeml_results.txt is a file produced from codeml under pairwise mode
Last modified December 12, 2020 
Jessie Pelosi

'''
import re
import sys 
from itertools import islice
args = sys.argv 
input_txt = sys.argv[1]

results = []
with open(input_txt, 'r') as infile:
	for line in infile:
		if "BMJR" in line and "NDUV" in line:
			results.append(list(islice(infile,4))[-1])
		output = ''.join(results)
		if "BMJR" in line and "SKYV" in line:
			results.append(list(islice(infile,4))[-1])
		output = ''.join(results)
		with open(input_txt + "ratios", "w") as outfile:
			outfile.write(output); outfile.close()

import csv
ratios = []
with open(input_txt + "ratios", 'r') as ratios_file:
		rd = csv.reader(ratios_file, delimiter = ' ')
		for row in rd:
			if "dN/dS=" in row[12]:
				ratios.append(row[14])
			elif "dN/dS=" in row[13]:
				ratios.append(row[15])
			elif "dN/dS=" in row[14]:
				ratios.append(row[16])
			elif "dN/dS=" in row[15]:
				ratios.append(row[17])
			elif "dN/dS=" in row[16]:
				ratios.append(row[18])
			elif "dN/dS=" in row[17]:
				ratios.append (row[19])
V_app_A_rad = ratios[0]
V_lin_A_rad = ratios[1]

with open("dn_ds_ratios.txt", "a") as outfile:
	outfile.write('\n'+input_txt+'\t'+str(V_app_A_rad)+'\t'+ str(V_lin_A_rad)); outfile.close()


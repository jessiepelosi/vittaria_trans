'''
Purpose: calculate GC content for each aligned cds, calculate the difference in GC content between aligned cds
Usage: python GC_content.py input_fasta 
Arguments:
	Input_fasta: aligned fasta file that you want to extract GC content from, contains at least one cds

Last Modidied 3 December 2022

** This should be used with the Vittaria-only ortholog set where there are only two sequences in each orthogroup **

'''

import sys 
args = sys.argv 
input_fasta = sys.argv[1]

input_file = open(input_fasta, 'r')
GC_freq_file = open(input_fasta+"gc.tsv", 'w')

from Bio import SeqIO
from Bio.SeqUtils import GC123 #function will write out %GC for total, first, second, and third codon positions of a sequence 

for cur_record in SeqIO.parse(input_file, "fasta"):
	contig_name = cur_record.name
	sequence = cur_record.seq
	content = str(GC123(sequence))
	con_temp1 = content.replace('(','') #function returns tuple, remove parentheses 
	con_temp2 = con_temp1.replace(')','')
	output_line = '%s\n' %con_temp2
	GC_freq_file.write(output_line)
input_file.close()
GC_freq_file.close()

import csv
GC_total=[]
GC_third=[]
with open(input_fasta+"gc.tsv", "rt") as tsv_file:
	rd = csv.reader(tsv_file)
	for row in rd:
		GC_total.append(row[0])
		GC_third.append(row[3]) 

GC_app_t = float(GC_total[0]) #V. appalachiana is always the first sequence in these ortholog pairs 
GC_lin_t = float(GC_total[1])
GC_app_3 = float(GC_third[0])
GC_lin_3 = float(GC_third[1])
difference_t = GC_app_t - GC_lin_t #calculate difference between V. appalachiana and V. lineata total GC content for each ortholog
difference_3 = GC_app_3 - GC_lin_3 #calculate difference between V. appalachiana and V. lineata 3rd codon position GC content for each ortholog 

with open("GC_dif.txt", "a") as outfile:
	outfile.write('\n'+input_fasta+'\t'+str(GC_app_t)+'\t'+ str(GC_lin_t)+'\t'+str(difference_t)+'\t'+str(GC_app_3)+'\t'+str(GC_lin_3)+'\t'+str(difference_3)); outfile.close()

outfile.close()

''' 
Purpose: run CodeMl for given alignments 
	Usage: python codeml.py input_fasta species_tree working_directory
	Arguments:
		Input_fasta: input gene alignment  
		Species_tree: input species tree 
		Working_directory: input the working directory for codeml to work from 
Script modified from https://biopython.org/wiki/PAML

Last modified 18 May 2020 
'''

import sys
args = sys.argv
input_fasta = sys.argv[1]
species_tree = sys.argv[2]
wk_dir = sys.argv[3]

input_file = open(input_fasta, "rt")

from Bio.Phylo.PAML import codeml
cml = codeml.Codeml()
cml.alignment = input_fasta
cml.tree = species_tree
cml.out_file = input_fasta + "_pairwise.txt"
cml.working_dir = wk_dir

cml.set_options(noisy = 9)
cml.set_options(verbose = 1)
cml.set_options(runmode = -2)
cml.set_options(seqtype = 1)
cml.set_options(CodonFreq = 2)
cml.set_options(ndata = 1)
cml.set_options(clock = 0)
cml.set_options(aaDist = 0)
cml.set_options(model = 1)
cml.set_options(NSsites = [0])
cml.set_options(icode = 0)
cml.set_options(Mgene = 0)
cml.set_options(fix_kappa = 0)
cml.set_options(kappa = 0)
cml.set_options(fix_omega = 0)
cml.set_options(omega = 0.4)
cml.set_options(fix_alpha = 1)
cml.set_options(alpha = 0)
cml.set_options(Malpha = 0)
cml.set_options(ncatG = 8)
cml.set_options(getSE = 0)
cml.set_options(RateAncestor = 0)
cml.set_options(Small_Diff = 5e-7)
cml.set_options(cleandata = 1)
cml.set_options(fix_blength = 0)
cml.set_options(method = 0)

cml.ctl_file = input_fasta + "control"
cml.write_ctl_file()

cml.run(ctl_file = input_fasta + "control", verbose = True)



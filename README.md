# <i> Vittaria appalachiana </i> comparative transcriptomics 
Code for "Life without a sporophyte: the population genetics and comparative transcriptomics of a gametophyte-only fern"- Pelosi, J.A., E.B. Sessa. 

# Orthology Determination and Alignment

Using OrthoFinder v. 2.3.11 (Emms and Kelly 2015, 2019) generate orthogroups as follows: 

<b> Set 1: </b> <i> Vittaria appalachiana </i> and <i> V. lineata </i> 

<b> Set 2: </b> <i> Adiantum raddianum, V. appalachiana, V. lineata </i>

<b> Set 3: </b> Pteridaceae + <i> Asplenium platyneuon </i> as an outgroup 

<TT> orthofinder -M msa -f proteomes/ </TT>

Retrieve single copy orthologs for OrthoFinder output and extract corresponding CDS (modify to use loop or array):

<TT> python extract_cds.py OG#######.pep transcriptome_1.fasta transcriptome_2.fasta transcriptome_n.fasta </TT>

Align CDS files using MACSE v. 2.04 (Ranwez et al. 2011) (modify to use loop or array):

<TT> java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq OG#######.cds </TT>

# GC Content Differences

Calculate GC content and difference between <i> V. appalachiana </i> and <i> V. lineata </i> across each ortholog and in the third codon position (modify to use loop or array). 

<TT> python GC_content.py OG#######.cds_NT </TT>

Load output file "GC_dif.txt" into R script "GC_diff.R" for comparison. 

# Calculate dN/dS ratios 

<b> Pairwise dN/dS: </b>  


Load output file "dn_ds_ratios.txt" into R script "dn_ds_diff.R" for comparion. 

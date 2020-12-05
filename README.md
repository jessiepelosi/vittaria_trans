# <i> Vittaria appalachiana </i> comparative transcriptomics 
Code for "Life without a sporophyte: the population genetics and comparative transcriptomics of a gametophyte-only fern"- Pelosi, J.A., E.B. Sessa. 

# Orthology Determination and Alignment

Using OrthoFinder v. 2.3.11 (Emms and Kelly 2015, 2019). 

<b> Set 1: </b> <i> Vittaria appalachiana </i> and <i> V. lineata </i> 

<b> Set 2: </b> <i> Adiantum raddianum, V. appalachiana, V. lineata </i>

<b> Set 3: </b> Pteridaceae + <i> Asplenium platyneuon </i> as an outgroup 

Retrieve single copy orthologs for OrthoFinder output and extract corresponding CDS:

<TT> python extract_cds.py OG#######.pep transcriptome_1.fasta transcriptome_2.fasta transcriptome_n.fasta </TT>

# GC Content Differences

<TT> for file in *.fasta; do python GC_content.py "$file"; done </TT> 

Load output file "GC_dif.txt" into R script "GC_diff.R" for comparison. 

# Calculate dN/dS ratios 

<b> Pairwise dN/dS: </b>  


Load output file "dn_ds_ratios.txt" into R script "dn_ds_diff.R" for comparion. 

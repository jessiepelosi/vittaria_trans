# <i> Vittaria appalachiana </i> comparative transcriptomics 
Code for "Life without a sporophyte: the population genetics and comparative transcriptomics of a gametophyte-only fern"- Pelosi, J.A., E.B. Sessa. 

All transcriptome assemblies are from the One Thousand Plant Transcriptome Initiatve (2019) and available at: https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019

# Orthology Determination and Alignment

Generate pep and CDS files using TransDecoder v. (Haas et al. 2014):

<TT> TransDecoder.LongOrfs -t [file.fasta] 
  TransDecoder.Predict -t [file.fasta] </TT>

Using OrthoFinder v. 2.3.11 (Emms and Kelly 2015, 2019) generate orthogroups as follows: 

<b> Set 1: </b> <i> Vittaria appalachiana </i> and <i> V. lineata </i> 

<b> Set 2: </b> <i> Adiantum raddianum, V. appalachiana, V. lineata </i>

<b> Set 3: </b> Pteridaceae + <i> Asplenium platyneuon </i> as an outgroup 

<TT> orthofinder -M msa -f proteomes/ </TT>

Retrieve single copy orthologs for OrthoFinder output and extract corresponding CDS (modify to use loop or array):

<TT> python extract_cds.py OG#######.pep transcriptome_1.cds transcriptome_2.cds transcriptome_n.cds </TT>

Align CDS files using MACSE v. 2.04 (Ranwez et al. 2011) (modify to use loop or array):

<TT> java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq OG#######.cds </TT>

Generate gene trees with IQTREE2 v. 2.1.0 (Minh et al. 2020, Hoang et al. 2018): 

<TT> iqtree2 -s OG#######.cds_NT --alrt 1000 -B 1000 -m TEST --redo </TT> 

Generate species tree with ASTRAL-III v. 5.6.2 (Zhang et al. 2018):

<TT> cat *.treefile > OGs_treefiles.tre
  
 astral -i OGs_treefiles.tre -o OGs_sp.tre 2> OGs_sp_tre.log </TT>  

# GC Content Differences

Calculate GC content and difference between <i> V. appalachiana </i> and <i> V. lineata </i> across each ortholog and in the third codon position (modify to use loop or array). 

<TT> python GC_content.py OG#######.cds_NT </TT>

Load output file "GC_dif.txt" into R script "GC_diff.R" for comparison. 

# dN/dS ratios 

Calculate pairwise dN/dS ratios using codeml in PAML v. 4.9h (Yang 2007) (modify to use loop or array).

<TT> python codeml.py OG#######.cds_NT species_tree.tre . </TT>

For this step, we used only three taxa, with the species tree: (Adiantum_radiannum,(Vittaria_appalachiana, Vittaria_lineata)); 

Next, extract the appropriate dN/dS ratios from the codeml results: 

<TT> python extract_paml_results_1.py OG#######.cds_NT_pairwise.txt </TT> 

Load output file "dn_ds_ratios.txt" into R script "dn_ds_diff.R" for comparison. 

# Gene Family Evolution

Load Orthogroups.GeneCount.tsv output file from OrthoFinder (for Pteridacaeae + <i> Asplenium </i>) into Count v. 9.1106 (Csűrös 2010). Run Count under Wagner Parsimony with  gain penatly of 1.2 (<i> Amborella </i> Genome Project 2013, Li et al. 2018) to identify gene families (orthogroups) that have expanded/contracted or gained/lost in <i> V. appalachiana </i>. Export Wagner Parsimony analysis from Count and Orthogroups.txt from OrthoFinder and input these into "gene_familiy_evol.R". 

# Inferring Whole Genome Duplications

<b> Ks plots </b> 

Use wgd (Zwaenepoel and Van de Peer 2019) to generate paralog age distributions (Ks) plots for <i> V. appalachiana </i> and <i> V. lineata </i>, the two known polyploids (based on chromosome counts, Gastony 1977). 

<TT> wgd mcl -s [NDUV/SKYV].cds --cds --mcl
    wgd ksd [NDUV/SKYV].mcl [NDUV/SKYV].cds
    wgd mix --method bgmm [NDUV/SKYV].cds.ks.tsv </TT>
  
Import the .tsv files into "ks_plots.R". 

<b> MAPS </b> 



<b> Determine topology of duplicates </b> 

Get orthogroups where there are exactly two copies of <i> V. appalachiana </i> and <i> V. lineata </i> and one copy of the outgroup <i> Asplenium platyneuron </i>. 

<TT> bash get_duplicate_trees.sh  </TT>

End up with 540 OGs. 

Extract CDS, align with MACSE, generate gene trees with IQTREE: 

<TT> python extract_cds.py OG#######.pep BMJR.cds FLTD.cds KJZG.cds POPJ.cds UJTT.cds WQML.cds YCKE.cds DCDT.cds GSXD.cds NDUV.cds SKYV.cds WCLG.cds XDDT.cds </TT> 
 java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq OG#######.cds </TT>
 iqtree2 -s OG#######.cds_NT --alrt 1000 -B 1000 -m TEST --redo </TT> 

Re-root trees with Newick Utilities v. 1.6 (Junier and Zdobnov 2010):

<TT> for file in *.treefile; do nw_reroot “$file” KJZG > “$file”_rerooted.tre; done </TT> 

Get orthogroups with certain topologies. 

<TT> grep -r "scaffold\-SKYV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_lin_v_app_dups.txt  #139, remove duplicates, ends up at 131 trees 
     grep -r "scaffold\-NUDV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_app_v_app_dups.txt #56, no duplicates 
     grep -r "scaffold\-NUDV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_app_v_lin_dups.txt #685, remove duplicates, ends up at 427 trees </TT>
 
Remove overlapping orthogroups in "v_lin_v_app_dups.txt" and "v_app_v_lin_dups.txt". In total, 457 trees support (((<i> V_appalachiana, V_lineata </i>), (<i> V_appalachiana, V_lineta </i>)), <i> Asplenium_platyneuron </i>), 56 support (((<i>V_appalachiana, V_appalachiana</i>), (<i>V_lineta, V_lineata</i>)),<i>Asplenium_platyneuron</i>). 

# <i> Vittaria appalachiana </i> comparative transcriptomics 
Code for  "Life without a sporophyte: the origin and consequences of asexual reproduction in a gametophyte-only fern" by Jessie A. Pelosi, W. Brad Barbazuk, and Emily B. Sessa. 

All transcriptome assemblies are from the One Thousand Plant Transcriptome Initiatve (2019) and available at: https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019

# Orthology Determination and Alignment

Generate pep and CDS files using TransDecoder v. 5.5 (Haas et al. 2014):
```
TransDecoder.LongOrfs -t [file.fasta] 

TransDecoder.Predict -t [file.fasta] 
```
Using OrthoFinder v. 2.3.11 (Emms and Kelly 2015, 2019) generate orthogroups as follows: 

<b> Set 1: </b> <i> Vittaria appalachiana </i> and <i> V. lineata </i> 

<b> Set 2: </b> <i> Adiantum raddianum, V. appalachiana, V. lineata </i>

<b> Set 3: </b> Pteridaceae + <i> Asplenium platyneuon </i> as an outgroup 
```
orthofinder -M msa -f proteomes/ 
```
Retrieve single copy orthologs for OrthoFinder output and extract corresponding CDS (modify to use loop or array):

```
python extract_cds.py OG#######.pep transcriptome_1.cds transcriptome_2.cds transcriptome_n.cds
```

Align CDS files using MACSE v. 2.04 (Ranwez et al. 2011) (modify to use loop or array):

```
java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq OG#######.cds 
```

Generate gene trees with IQTREE2 v. 2.1.0 (Minh et al. 2020, Hoang et al. 2018):

```
iqtree2 -s OG#######.cds_NT --alrt 1000 -B 1000 -m TEST --redo 
```

Generate species tree with ASTRAL-III v. 5.6.2 (Zhang et al. 2018):
```
cat *.treefile > OGs_treefiles.tre 
  
astral -i OGs_treefiles.tre -o OGs_sp.tre 2> OGs_sp_tre.log 
```

# GC Content Differences

Calculate GC content and difference between <i> V. appalachiana </i> and <i> V. lineata </i> across each ortholog and in the third codon position (modify to use loop or array). 

```
python GC_content.py OG#######.cds_NT
```

Load output file "GC_dif.txt" into R script "GC_diff.R" for comparison. 

# dN/dS ratios 

Calculate pairwise dN/dS ratios using codeml in PAML v. 4.9h (Yang 2007) (modify to use loop or array).

```
python codeml.py OG#######.cds_NT species_tree.tre .
```

For this step, we used only three taxa, with the species tree: (Adiantum_radiannum,(Vittaria_appalachiana, Vittaria_lineata)); 

Next, extract the appropriate dN/dS ratios from the codeml results: 

```
python extract_paml_results_1.py OG#######.cds_NT_pairwise.txt 
```

Load output file "dn_ds_ratios.txt" into R script "dn_ds_diff.R" for comparison. 

# Gene Family Evolution

Load Orthogroups.GeneCount.tsv output file from OrthoFinder (for Pteridacaeae + <i> Asplenium </i>) into Count v. 9.1106 (Csűrös 2010). Run Count under Wagner Parsimony with  gain penatly of 1.2 (<i> Amborella </i> Genome Project 2013, Li et al. 2018) to identify gene families (orthogroups) that have expanded/contracted or gained/lost in <i> V. appalachiana </i>. Export Wagner Parsimony analysis from Count and Orthogroups.txt from OrthoFinder and input these into "gene_familiy_evol.R". 

# Inferring Whole Genome Duplications

<b> Ks plots </b> 

Use wgd (Zwaenepoel and Van de Peer 2019) to generate paralog age distributions (Ks) plots for <i> V. appalachiana </i> and <i> V. lineata </i>, the two known polyploids (based on chromosome counts, Gastony 1977). 

```
wgd mcl -s [NDUV/SKYV].cds --cds --mcl 

wgd ksd [NDUV/SKYV].mcl [NDUV/SKYV].cds

wgd mix --method bgmm [NDUV/SKYV].cds.ks.tsv
```
  
Import the .tsv files into "ks_plots.R". 

<b> MAPS </b> 

Select taxa for MAPS analysis (tree must be ladderized, thus several sister taxa must be removed from the main species tree). Then select orthogroups where at least one gene copy is present for each taxon. MAPS is described in Li et al. 2015 and Li et al. 2018, and the source code can be found here: https://bitbucket.org/barkerlab/maps/src/master/

```
bash get_MAPS_OGs.sh
```

Create gene trees for each orthogroup.  

```
python extract_cds.py OG#######.pep FLTD.cds KJZG.cds WQML.cds DCDT.cds NDUV.cds SKYV.cds WCLG.cds  
 
java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq OG#######.cds

iqtree2 -s OG#######.cds_NT --alrt 1000 -B 1000 -m TEST --redo 
```

Rename tip labels, remove bootstrap support values, and concatenate into one file (pteridaceae_noBS.tre). 

Run MAPS with the list file (ladderized species tree) and gene tree file. 

```
perl maps.pl --l pteridaceae.list --t pteridaceae_noBS.tre
```
<b> For MAPS Simulations: </b>

Calculate geometric mean at root and estimate birth and death rates of gene families in `MAPS.R`- this will prepare the input for running null and positive simulations with MAPS. Some scripts have been slightly modified from Li et al. 2015, 2018 to change the version of JPRIME being called. 

```
perl simulateGeneTrees.3.0_JAP.pl #change sim.ctl for each set of birth and death rates 
```
For null simulations: 
```
perl sampleTrees.pl -in NULL3000trees.txt -n 1000 -r 100 -out NULL_subsamp

for file in NULL_subsamp.*.tre; do perl maps.pl --l pteridaceae.list --tree "$file" --o NULL4000_"$file"_out;done
```
For positive simulations: 
```
perl sampleTrees.pl -in POS3000trees.txt -n 1000 -r 100 -out POS_subsamp

for file in POS_subsamp.*.tree; do perl maps.pl --l pteridaceae.list --tree "$file" --o POS4000_"$file"_out;done
```

Aggregate data with awk for each node (change number of commands depending on number of nodes in analysis) and run Fisher's exact tests. 

```
awk -F , 'FNR == 2{print $3}' 4000_subsamp.*/*.csv >> MAPS_simulations.txt
awk -F , 'FNR == 3{print $3}' 4000_subsamp.*/*.csv >> MAPS_simulations.txt
awk -F , 'FNR == 4{print $3}' 4000_subsamp.*/*.csv >> MAPS_simulations.txt
awk -F , 'FNR == 5{print $3}' 4000_subsamp.*/*.csv >> MAPS_simulations.txt
awk -F , 'FNR == 6{print $3}' 4000_subsamp.*/*.csv >> MAPS_simulations.txt
```

Summarize MAPS analyses on simulated data and then run Fisher's exact tests. 

```
cd null_sims/

python summarize_MAPS.py 5 1000 Pteridaceae

perl runFisher_null.pl

cd ../pos_sims/

python summarize_MAPS.py 5 1000 Pteridaceae 

perl runFisher_positive.pl 
```

<b> Determine topology of duplicates </b> 

Get orthogroups where there are exactly two copies of <i> V. appalachiana </i> and <i> V. lineata </i> and one copy of the outgroup <i> Asplenium platyneuron </i>. Given the placement of WGDs suggested by MAPs, we also examined the topology of duplicates where there were two copies of <i> V. appalachiana </i> and <i> V. lineata </i> and one copy of <i> Adiantum aleuticum </i>. 

<TT> bash get_duplicate_trees.sh  #or bash get_duplicate_trees_4.sh for Adiantum</TT>

End up with 540 OGs with <i> Asplenium </i> as the outgroup, and 471 OGs with <i> Adiantum </i> as the outgroup. 

Extract CDS, align with MACSE, generate gene trees with IQTREE: 
```
python extract_cds.py OG#######.pep BMJR.cds FLTD.cds KJZG.cds POPJ.cds UJTT.cds WQML.cds YCKE.cds DCDT.cds GSXD.cds NDUV.cds SKYV.cds WCLG.cds XDDT.cds 
 
java -Xms4000m -jar macse_v2.04.jar -prog alignSequences -seq OG#######.cds

iqtree2 -s OG#######.cds_NT --alrt 1000 -B 1000 -m TEST --redo 
```
Re-root trees with Newick Utilities v. 1.6 (Junier and Zdobnov 2010):
```
for file in *.treefile; do nw_reroot “$file” KJZG > “$file”_rerooted.tre; done #Use WCLG for Adiantum as outgroup 
```
Get orthogroups with certain topologies. 

For <i> Asplenium </i> as outgroup: 

```
grep -r "scaffold\-SKYV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_lin_v_app_dups.txt  #139, remove duplicates, ends up at 131 trees 

grep -r "scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_app_v_app_dups.txt #56, no duplicates

grep -r "scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-SKYV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_app_v_lin_dups.txt #685, remove duplicates, ends up at 427 trees 
``` 
Remove overlapping orthogroups in "v_lin_v_app_dups.txt" and "v_app_v_lin_dups.txt". In total, 457 trees support (((<i> V. appalachiana, V. lineata </i>), (<i> V. appalachiana, V. lineta </i>)), <i> Asplenium platyneuron </i>), 56 support (((<i>V. appalachiana, V. appalachiana</i>), (<i>V. lineta, V. lineata</i>)),<i>Asplenium platyneuron</i>). 

For <i> Adiantum </i> as outgroup: 

``` 
grep -r "scaffold\-SKYV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_lin_v_app_dups.txt  #130, remove duplicates, ends up at 120 trees

grep -r "scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_app_v_app_dups.txt #53, no duplicates 

grep -r "scaffold\-NDUV\-[A-Za-z0-9\.\_\:\-]*,scaffold\-SKYV\-[A-Za-z0-9\.\_\:\-]*" *.tre -o > v_app_v_lin_dups.txt #596, remove duplicates, ends up at 367 trees 
``` 
Remove overlapping orthogroups in "v_lin_v_app_dups.txt" and "v_app_v_lin_dups.txt". In total, 396 trees support (((<i> V. appalachiana, V. lineata </i>), (<i> V. appalachiana, V. lineta </i>)), <i> Adiantum aleuticum</i>), 53 support (((<i>V. appalachiana, V. appalachiana</i>), (<i>V. lineta, V. lineata</i>)),<i> Adiantum aleuticum </i>). 

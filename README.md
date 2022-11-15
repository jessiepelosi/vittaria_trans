# <i> Vittaria appalachiana </i> comparative transcriptomics 
Code for  "Life without a sporophyte: the origin and consequences of asexual reproduction in a gametophyte-only fern" by Jessie A. Pelosi, W. Brad Barbazuk, and Emily B. Sessa. 

All transcriptome assemblies are from the One Thousand Plant Transcriptome Initiatve (2019) and available at: https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019

This repository is divided into four folders: 
1. [GC_content](https://github.com/jessiepelosi/vittaria_trans/tree/main/GC_Content)
2. [dN/dS](https://github.com/jessiepelosi/vittaria_trans/tree/main/dNdS)
3. [genefams](https://github.com/jessiepelosi/vittaria_trans/tree/main/genefams)
4. [WGD](https://github.com/jessiepelosi/vittaria_trans/tree/main/WGD)

The process for generating of orthologs for use in downstream anlayses is below. 

# Orthology Determination and Alignment

Generate pep and CDS files using TransDecoder v. 5.5 (Haas et al. 2014):
```
TransDecoder.LongOrfs -t [file.fasta] 

TransDecoder.Predict -t [file.fasta] 
```
Using OrthoFinder v. 2.3.11 (Emms and Kelly 2015, 2019) generate orthogroups as follows: 

<b> Set 1: </b> <i> Adiantum raddianum, V. appalachiana, V. lineata </i>

<b> Set 2: </b> Pteridaceae + <i> Asplenium platyneuon </i> as an outgroup 
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

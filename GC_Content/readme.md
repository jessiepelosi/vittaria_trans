# GC Content Differences

Calculate observed GC content and difference between <i> V. appalachiana </i> and <i> V. lineata </i> across each ortholog and in the third codon position (modify to use loop or array). 

```
python GC_content.py OG#######.cds_NT
```

Load output file "GC_dif.txt" into R script "GC_diff.R" for comparison. 

Equilibrium GC content (GC* and GC3*) were calculated with BppML in the Bio++Suite ver 2.3.1. An initial tree was generated with IQTREE with the single-copy orthologs between <i>Vittaria appalachiana, V. lineata,</i> and <i>Adiantum raddianum</i> as a starting point for branch length estimation. The non-homogeneous model proposed by Galtier and Gouy (1998) was used as a modification of the Tamaura (1992) model of sequence evolution, allowing each branch to have a different equilibirium GC content and a fixed transition/transversion ratio. A gamma distribution with 4 rate classes and invariant site rate distribution was used. 
```
bppml param="$file".ctl
```
<b>Include config file.</b>

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

# Reptitive Elements

We used two approaches (assembly and read-based) to assess repeat content in the <i>Vittaria appalachiana</i> and <i>V. lineata</i> transcriptomes. 

## RepeatMasker 

RepeatMasker ver. 4.0.5 was used with a Viridiplantae repeat set for the assembly-based approach. 

```

```

## dnaPipeTE

dnaPipeTE ver 1.4c was used as the reads-based approach. The singularity container was downloaded from [here](https://github.com/clemgoub/dnaPipeTE) and run in a singularity instance on the forward reads of each library (downloaded from 1KP). 

```
singularity shell --bind [project]:/mnt [path/to/dnaPipeTE.img]
# Once in singularity 
cd /opt/dnaPipeTE/
python3 dnaPipeTE.py -input /mnt/ERR204094[0/1]_1.fastq -output /mnt/output/ -RM_lib /mnt/RepeatMasker.lib -sample_size [no. reads] -RM_t 0.1 -cpu 6
```

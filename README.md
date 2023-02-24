## Project ticket number: 41

### Estimation of transposition frequencies

#### To estimate transposition frequency from minion reads for sample name "BC1" with a fastq file "BC1.fastq, run the following command in biowulf:
```bash
# Add sample name
SAMPLE=BC1

# Add path to directory containing fastq file
FASTQDIR=/my/path/to/reads/dir

# Add info about the type of vector, either 'original' or 'modified'
TYPE=original

# Run the pipeline like so:
run_IS_mapper.sh ${SAMPLE} ${FASTQDIR} ${TYPE}
```

#### This will generate the following files witin the 03-blastn directory:
```
BC1_trimmed.filter.hits.bed
BC1_trimmed.filter.hits.flanking_vs_original.NR.bn
```
Note that the name of the files will change depending on the sample name and the type of plasmid specified (original or modified).

#### Import the two files into the python jupyter notebook *predict_IS_from_single_end_blast.ipynb* to plot frequency of transposition per TTAA site

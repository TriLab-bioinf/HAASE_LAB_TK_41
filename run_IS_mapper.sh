#!/usr/bin/bash
set -o errexit

# Get sample prefix
SAMPLE=$1
FASTQDIR=$2
TYPE=$3 # 'original' or 'modified'

# Piggybac transposon prefix
PIGGYBAC=pTet_pBac_LE35_LE35

# Load required modules
module load seqtk mummer blast porechop

# OPTIONAL: Demultiplexing with readuck:
#readucks -t 14 -i ${SAMPLE}.fastq -o demuxed -b --pcr_barcodes --verbosity 1 --check_reads 10000 > readucks.log 2>&1 &

#Remove the minion adapters and bean reads based on barcodes with porechop like so:

## To trim:
TRIMDIR='01-trim'
if [[ ! -d ${TRIMDIR} ]]; then
    mkdir ${TRIMDIR}
fi

echo STEP-1: Trimming with porechop
porechop -i ${FASTQDIR}/${SAMPLE}.fastq -o ./${TRIMDIR}/${SAMPLE}_trimmed.fastq -t 14 > porechop.log 2>&1 &

## To demultiplex:
# porechop -i ${SAMPLE}.fastq -t 14 -b 01-trimmed > porechop.log 2>&1

# Convert trimmed fastq sequences to fasta with seqtk:
seqtk seq -A ./${TRIMDIR}/${SAMPLE}_trimmed.fastq > ./${TRIMDIR}/${SAMPLE}_trimmed.fasta


# Align fasta reads to PiggyBac with nucmer:
echo STEP-2: Aligning fasta reads to PiggyBac with nucmer
MUMMERDIR='02-mummer'
if [[ ! -d ${MUMMERDIR} ]]; then
    mkdir ${MUMMERDIR}
fi
nucmer --mum -c 10 --minalign=80 --minmatch=10 -p ./${MUMMERDIR}/${SAMPLE}_vs_${PIGGYBAC} \
    ${PIGGYBAC}.fa ./${TRIMDIR}/${SAMPLE}_trimmed.fasta

# Generate coordinate file from delta alignment:
show-coords -clrTH ./${MUMMERDIR}/${SAMPLE}_vs_${PIGGYBAC}.delta \
    > ./${MUMMERDIR}/${SAMPLE}_vs_${PIGGYBAC}_filter.coords

# Extract read IDs matching Piggybac from coords file:
cat ./${MUMMERDIR}/${SAMPLE}_vs_${PIGGYBAC}_filter.coords | \
perl -e 'while(<>){chomp;@x=split /\t/;$h{$x[12]}+=$x[9] };foreach $k (keys %h){print "$k\n" if $h{$k}>80 }' \
> ./${MUMMERDIR}/${SAMPLE}_vs_${PIGGYBAC}_hits.ids

# OPTIONAL: Generate histogram of read-length distribution from coords file for those reads matching PiggyBac:
#grep -wf ${SAMPLE}_vs_${PIGGYBAC}_hits.ids ${SAMPLE}_vs_${PIGGYBAC}_filter.coords | \
#cut -f 9 | binner.pl 100 >& ${SAMPLE}_vs_${PIGGYBAC}_filter.hits.length.hist

# Extract fasta reads matching Piggybac sequence from list of IDs:
seqtk subseq ./${TRIMDIR}/${SAMPLE}_trimmed.fasta ./${MUMMERDIR}/${SAMPLE}_vs_${PIGGYBAC}_hits.ids > ./${MUMMERDIR}/${SAMPLE}_trimmed.hits.fasta

# Make blastn database with fasta file of reads that hit Piggybac
makeblastdb -in ./${MUMMERDIR}/${SAMPLE}_trimmed.hits.fasta -dbtype nucl

# Alternative code for extracting single-end flanking seqs
echo STEP-3: Run blastn to extract flanking sequences
BLASTNDIR='03-blastn'
if [[ ! -d ${BLASTNDIR} ]]; then
    mkdir ${BLASTNDIR}
fi

blastn -query ${PIGGYBAC}.fa -db ./${MUMMERDIR}/${SAMPLE}_trimmed.hits.fasta -evalue 1e-10 \
    -outfmt '7 qaccver saccver pident length mismatch gapopen qstart qend sstart send qlen slen evalue' \
    -max_target_seqs 100000 -out ./${BLASTNDIR}/${SAMPLE}_trimmed.hits.bn

./filter_blastn.pl < ./${BLASTNDIR}/${SAMPLE}_trimmed.hits.bn > ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.bed

# Extract flanking sequences using bed file from the reads:
seqtk subseq ./${MUMMERDIR}/${SAMPLE}_trimmed.hits.fasta ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.bed > ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking.fasta

# OPTIONAL: Remove duplicated fasta sequences (if there is any):
# cat ${SAMPLE}_trimmed.filter.hits.flanking.fasta|perl -ne 'if(m/^(>\S+)/){$h{$1}++}else{print "$1\n$_" if $h{$1}==1}' > ${SAMPLE}_trimmed.filter.hits.flanking.dedup.fasta

# Blast flanking sequences against plasmid to extract exact coordinates:
echo STEP-4:  Blastn flanking sequences against plasmid to extract coordinates of IS

## Generate blastn DBs
if [[ ${TYPE} == 'original' ]]; then
    if [[ ! -f pHSG298_original_FIXED.fa.njs ]]; then
        makeblastdb -in pHSG298_original_FIXED.fa -dbtype nucl
    fi    
    ## For original plasmid:
    blastn -query ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking.fasta \
        -db pHSG298_original_FIXED.fa -max_target_seqs 1 -evalue 1e-10 -word_size 11 \
        -gapopen 5 -gapextend 2 -penalty -3 -reward 2 \
        -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send qlen slen evalue' \
        > ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking_vs_original.bn

    # Remove duplicated hits (keep the one with the longest length)
    cat  ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking_vs_original.bn|perl -e 'while(<>){@x=split /\t/; if ($x[3] > $h{$x[0]}){$L{$x[0]} = $_; $h{$x[0]} = $x[3]}  }foreach $k (keys %L){print "$L{$k}"}' \
        > ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking_vs_original.NR.bn

else
    if [[ ! -f pHSG298plus8NdeI_modified_FIXED.fa.njs ]]; then
        makeblastdb -in pHSG298plus8NdeI_modified_FIXED.fa -dbtype nucl
    fi
    ## For modified plasmid:
    blastn -query ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking.fasta \
        -db pHSG298plus8NdeI_modified_FIXED.fa -max_target_seqs 1 -evalue 1e-10 -word_size 11 \
        -gapopen 5 -gapextend 2 -penalty -3 -reward 2 \
        -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send qlen slen evalue' \
        > ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking_vs_modified.bn

    cat  ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking_vs_modified.bn|perl -e 'while(<>){@x=split /\t/; if ($x[3] > $h{$x[0]}){$L{$x[0]} = $_; $h{$x[0]} = $x[3]}  }foreach $k (keys %L){print "$L{$k}"}' \
        > ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking_vs_modified.NR.bn

fi

# Upload ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.flanking_vs_original.NR.bn and 
# ./${BLASTNDIR}/${SAMPLE}_trimmed.filter.hits.bed into jupyter-lab for further analysis
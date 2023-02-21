#!/usr/bin/bash
set -o errexit

SAMPLE=$1
BLASTNDIR=$2

# Load required modules
module load seqtk clustalw cd-hit

# Generate bed file for each IS site +/- 50bp
cat ./${BLASTNDIR}/${SAMPLE}_trimmed.hits.bn | \
    perl -ne '@x=split /\t/;($x[8],$x[9])=($x[8],$x[9]) if $x[8]>$x[9];$s5=$x[8]-52;$e5=$x[8]+50;$e3=$x[9]+50; $s3=$x[9]-50;print "$x[1]\t$s5\t$e5\n$x[1]\t$s3\t$e3\n"' \
    > ${SAMPLE}_trimmed.filter.hits.PLPD.flank.bed

# Extract ~100bp sequences spanning either IS site of the read with seqtk
seqtk subseq ./${MUMMERDIR}/${SAMPLE}_trimmed.hits.fasta ${SAMPLE}_trimmed.filter.hits.PLPD.flank.bed \
    > ${SAMPLE}_trimmed.filter.hits.PLPD.flank.fasta

# Note-1: Edit fasta headers by removing long prefixes (e.g.'Hickman_n2n_3_1A71') to make headers short enough so they are not truncated in the output of cd-hit-est

# To generate Logo plots describing conservation of IS sequences:

# Cluster flanking sequences spanning either IS site of the read with cd-hit-est:
cd-hit-est -i ${SAMPLE}_trimmed.filter.hits.PLPD.flank.fasta -o ${SAMPLE}_trimmed.filter.hits.PLPD.flank.cdhit -c 0.9 -gap –3 -g 1

# Create output directory:
CDHITDIR=${SAMPLE}_CDHIT
if [[ ! -f ${CDHITDIR} ]]; then
    mkdir ${CDHITDIR}
fi

# Extract fasta sequences for all reads within each cluster:
./extract_cluster_sequences.pl -c ${SAMPLE}_trimmed.filter.hits.PLPD.flank.cdhit.clstr \
    -f ${SAMPLE}_trimmed.filter.hits.PLPD.flank.fasta \
    -p ./${CDHITDIR}/${SAMPLE}

# Go to sample directory:
cd ${CDHITDIR}

# Forward PB 5'
SEQ=ATTATCATGACATTAACCTAT

rm –f ${SAMPLE}_5p.fasta

for i in `grep ${SEQ} *_Cluster_*.fasta | cut -f 1 -d ':' | sort -u`; do cat $i \
    >> ${SAMPLE}_5p.fasta
done

# Reverse PB 5'
SEQ2=ATAGGTTAATGTCATGATAAT

for i in `grep ${SEQ2} *_Cluster_*.fasta | cut -f 1 -d ':' | sort -u`; do cat $i >> ${SAMPLE}_5p_rev.fasta; done

cat ${SAMPLE}_5p_rev.fasta | \
    perl -ne 'if(m/^>/){print}else{$_=~tr/ACGTN/TGCAN/;$seq=reverse($_);print "$seq\n"}' \
    >> ${SAMPLE}_5p_revcomp.fasta

cat ${SAMPLE}_5p_revcomp.fasta >> ${SAMPLE}_5p.fasta

# Forward PB 3'
SEQ3=TCGCTAACGGATTCACCAC

rm –f ${SAMPLE}_3p.fasta

for i in `grep ${SEQ3} *_Cluster_*.fasta | cut -f 1 -d ':' | sort -u`; do cat $i \
    > ${SAMPLE}_3p.fasta; done

# Reverse PB 3'
SEQ4=GTGGTGAATCCGTTAGCGA

for i in `grep ${SEQ4} *_Cluster_*.fasta | cut -f 1 -d ':' | sort -u`; do cat $i \
    >> ${SAMPLE}_3p_rev.fasta 
done

cat ${SAMPLE}_3p_rev.fasta | \
    perl -ne 'if(m/^>/){print}else{$_=~tr/ACGTN/TGCAN/;$seq=reverse($_);print "$seq\n"}' \
    >> ${SAMPLE}_3p_revcomp.fasta

cat ${SAMPLE}_3p_revcomp.fasta >> ${SAMPLE}_3p.fasta

# Generate multiple sequence alignment of flanking seqs:

clustalw -INFILE=${SAMPLE}_5p.fasta -ALIGN -TREE
clustalw -INFILE=${SAMPLE}_3p.fasta -ALIGN -TREE






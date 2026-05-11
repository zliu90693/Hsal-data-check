#!/bin/bash

echo 'usage: ./dropseq.sh <R1 fastq> <R2 fastq> <sample name>'


#input is two fastq files from the dropseq library you want to analyze
R1=$1
R2=$2
sample_name=$3
dropseq_code_dir=/home/emily/bin/Drop-seq_tools-2.0.0/
picard=/home/emily/bin/Drop-seq_tools-2.0.0/picard.jar
star_dir=/home/emily/data/genome_annotation/genomes/STAR_indexes/hsal_ncbi_v8.5_and_dmel
ref_fasta=/home/emily/data/genome_annotation/genomes/STAR_indexes/hsal_ncbi_v8.5_and_dmel/hsal_ncbi_v8.5_and_dmel.fasta
refflat=/home/emily/data/genome_annotation/annotation/hsal_ncbi_v8.5/hsal_ncbi_v8.5_krh1.full.refFlat

mkdir ${sample_name}_intermediate/
mkdir ${sample_name}_summaries/

java -jar $picard FastqToSam F1=$R1 F2=$R2 O=${sample_name}.bam SM=$sample_name

${dropseq_code_dir}/TagBamWithReadSequenceExtended INPUT=${sample_name}.bam OUTPUT=${sample_name}.unaligned_tagged_Cell.bam SUMMARY=${sample_name}.unaligned_tagged_Cellular.bam BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1

${dropseq_code_dir}/TagBamWithReadSequenceExtended INPUT=${sample_name}.unaligned_tagged_Cell.bam OUTPUT=${sample_name}.unaligned_tagged_CellMolecular.bam SUMMARY=${sample_name}.unaligned_tagged_Molecular.bam BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

mv ${sample_name}.bam ${sample_name}_intermediate/
mv ${sample_name}.unaligned_tagged_Cell.bam ${sample_name}_intermediate/

${dropseq_code_dir}/FilterBam TAG_REJECT=XQ INPUT=${sample_name}.unaligned_tagged_CellMolecular.bam OUTPUT=${sample_name}.unaligned_tagged_filtered.bam

mv ${sample_name}.unaligned_tagged_CellMolecular.bam ${sample_name}_intermediate/
mv ${sample_name}.unaligned_tagged_Cellular.bam ${sample_name}_intermediate/
mv ${sample_name}.unaligned_tagged_Molecular.bam ${sample_name}_intermediate/

${dropseq_code_dir}/TrimStartingSequence INPUT=${sample_name}.unaligned_tagged_filtered.bam OUTPUT=${sample_name}.unaligned_tagged_trimmed_smart.bam OUTPUT_SUMMARY=${sample_name}.adapter_trimming_report.txt SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5

mv ${sample_name}.unaligned_tagged_filtered.bam ${sample_name}_intermediate/

${dropseq_code_dir}/PolyATrimmer INPUT=${sample_name}.unaligned_tagged_trimmed_smart.bam OUTPUT=${sample_name}.unaligned_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=${sample_name}.polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6 USE_NEW_TRIMMER=true

mv ${sample_name}.unaligned_tagged_trimmed_smart.bam ${sample_name}_intermediate/

java -jar -Xmx8g $picard SamToFastq INPUT=${sample_name}.unaligned_mc_tagged_polyA_filtered.bam FASTQ=${sample_name}.unaligned_mc_tagged_polyA_filtered.fastq/

STAR --genomeDir $star_dir --readFilesIn ${sample_name}.unaligned_mc_tagged_polyA_filtered.fastq --outFileNamePrefix ${sample_name}

mkdir ${sample_name}_starLogs
mv ${sample_name}Log* ${sample_name}_starLogs
mv ${sample_name}SJ* ${sample_name}_starLogs

mv ${sample_name}.unaligned_mc_tagged_polyA_filtered.fastq ${sample_name}_intermediate/

java -Xmx8g -jar $picard SortSam I=${sample_name}Aligned.out.sam O=${sample_name}.aligned.sorted.bam SO=queryname

mv ${sample_name}Aligned.out.sam ${sample_name}_intermediate/

java -Xmx8g -jar $picard MergeBamAlignment REFERENCE_SEQUENCE=${ref_fasta} UNMAPPED_BAM=${sample_name}.unaligned_mc_tagged_polyA_filtered.bam ALIGNED_BAM=${sample_name}.aligned.sorted.bam OUTPUT=${sample_name}.merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

mv ${sample_name}.unaligned_mc_tagged_polyA_filtered.bam ${sample_name}_intermediate
mv ${sample_name}.aligned.sorted.bam ${sample_name}_intermediate/

${dropseq_code_dir}/TagReadWithGeneFunction I=${sample_name}.merged.bam O=${sample_name}.star_gene_exon_tagged.bam ANNOTATIONS_FILE=$refflat

mv ${sample_name}.merged.bam ${sample_name}_intermediate/

mv ${sample_name}*txt ${sample_name}_summaries
mv ${sample_name}_starLogs ${sample_name}_summaries


#!/bin/bash

WORKDIR="/mnt/e/workspace/bd_RNASeq/"
GENOME="/mnt/e/workspace/reference_genome/Bomo_genome_assembly.fa"
GTF="/mnt/e/workspace/reference_genome/Bomo_gene_models.gtf"
INDEX="/mnt/e/workspace/reference_genome/STAR_genome_index"
RSEMINDEX="/mnt/e/workspace/reference_genome/RSEM_index"
FILE_PROCESS_LIST="/mnt/e/workspace/bd_RNASeq/RNA_processing_list.txt"
DICT="/mnt/e/workspace/reference_genome/Bomo_genome_assembly.dict"


###indexing###
#This analysis was performed using Windows WSL on Windows partitions 
#In STAR, if run partition does not support FIFO (e.g. Windows partitions FAT, NTFS), it is required to point --outTmpDir to a Linux partition.
STAR --runThreadN 5 --runMode genomeGenerate --genomeDir $INDEX --genomeFastaFiles $GENOME --sjdbGTFfile $GTF --genomeSAindexNbases 13 --outTmpDir /home/tommy/STAR_temp
rsem-prepare-reference --gtf $GTF -p 4 $GENOME $RSEMINDEX
gatk CreateSequenceDictionary -R ${GENOME} -O ${DICT}
samtools faidx ${GENOME}

pushd $WORKDIR

for i in {1..3}; do
    READ_1=$(awk -v ind=$i 'NR==ind{print $1}' $FILE_PROCESS_LIST)
    READ_2=$(awk -v ind=$i 'NR==ind{print $2}' $FILE_PROCESS_LIST)
    PREFIX=$(awk -v ind=$i 'NR==ind{print $3}' $FILE_PROCESS_LIST)
    
    #READ_1 and READ_2 were already trimmed by Trim Galore! (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
    
    echo ${PREFIX}
    echo ${READ_1}
    echo ${READ_2}
    OPDIR=${PREFIX}_processed
    mkdir ${OPDIR}

    STAR --outFileNamePrefix ${OPDIR}/${PREFIX} --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outSAMstrandField intronMotif --sjdbGTFfile $GTF --outTmpDir /home/tommy/${PREFIX}_STAR_temp --outSAMattrRGline ID:$PREFIX CN:IGB LB:PairedEnd PL:Illumina PU:Unknown SM:$PREFIX --genomeDir ${INDEX} --runThreadN 5 --readFilesIn $READ_1 $READ_2 --twopassMode Basic
    rsem-calculate-expression --alignments --paired-end --estimate-rspd --strandedness reverse --no-bam-output -p 4 ${OPDIR}/${PREFIX}Aligned.toTranscriptome.out.bam $RSEMINDEX ${OPDIR}/${PREFIX}_rsem

    #sambamba sort
    sambamba sort -o ${OPDIR}/${PREFIX}_sorted.bam -p -t 5 ${OPDIR}/${PREFIX}Aligned.out.bam
    rm ${OPDIR}/${PREFIX}Aligned.out.bam
    #sambamba index
    sambamba index -p -t 5 ${OPDIR}/${PREFIX}_sorted.bam
    
    #mark duplicates
    gatk MarkDuplicates -I ${OPDIR}/${PREFIX}_sorted.bam -O ${OPDIR}/${PREFIX}_dupMarked.bam -M ${OPDIR}/${PREFIX}_dup.metrics -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT 2>${OPDIR}/${PREFIX}_MarkDuplicates.log

    #SplitNCigarReads
    gatk SplitNCigarReads -R ${GENOME} -I ${OPDIR}/${PREFIX}_dupMarked.bam -O ${OPDIR}/${PREFIX}_split.bam 2>${OPDIR}/${PREFIX}.SplitNCigarReads.log
    rm ${OPDIR}/${PREFIX}_dupMarked.bam
    rm ${OPDIR}/${PREFIX}_dupMarked.bai

    #Run HaplotypeCaller
    gatk HaplotypeCaller -R ${GENOME} -I ${OPDIR}/${PREFIX}_split.bam -O ${OPDIR}/${PREFIX}.vcf 2>${OPDIR}/${PREFIX}.HaplotypeCaller.log

    #Filter variants
    gatk VariantFiltration -R ${GENOME} -V ${OPDIR}/${PREFIX}.vcf -window 35 -cluster 3 --filter-name "FS" -filter "FS > 30.0" --filter-name "QD" -filter "QD < 2.0" -O ${OPDIR}/${PREFIX}_filtered.vcf 2>${OPDIR}/${PREFIX}.VariantFilter.log
done

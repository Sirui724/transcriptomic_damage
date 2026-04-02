#!/bin/bash 

module load SAMTools/1.10
module load bedtools/2.30.0
module load STAR/2.7.10a
module load rsem/1.3.3
module load perl
module load python/3.13.5


file_pathway=/data/gladyshev/sz1012/blood_data   #need change according to the fastq file path 
output_pathway=/data/gladyshev/sz1012/damage_project/out
star_reference_pathway=/data/gladyshev/sz1012/reference/STAR_index_GRCh38
star_reference_pathway2=/data/gladyshev/sz1012/reference/RSEM_index_GRCh38/RSEM_index_GRCh38
info_file=info_7.log

samples=(SRR9666181
)


for sample_name in ${samples[@]}
do
    echo "${sample_name}"

    ### 1st STAR run
    STAR --runThreadN 25 --runMode alignReads --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir ${star_reference_pathway} \
        --readFilesIn ${file_pathway}/${sample_name}_1.fastq,${file_pathway}/${sample_name}_2.fastq \
        --outFileNamePrefix ${output_pathway}/${sample_name} \
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 --chimOutType WithinBAM HardClip \
        --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 --chimMultimapNmax 50  

    ### 2nd STAR run
    STAR --runThreadN 25 \
        --genomeDir ${star_reference_pathway} \
        --readFilesIn ${file_pathway}/${sample_name}_1.fastq,${file_pathway}/${sample_name}_2.fastq \
        --outFileNamePrefix ${output_pathway}/${sample_name}_re \
        --outSAMmultNmax 1 --winAnchorMultimapNmax 200 \
        --alignIntronMax 100000 --alignMatesGapMax 100000 \
        --outSAMprimaryFlag OneBestScore --outSAMtype SAM

    ### Rename intermediate files
    samtools view -q 12 -b ${output_pathway}/${sample_name}Aligned.sortedByCoord.out.bam > ${output_pathway}/${sample_name}_uniqmap.bam
    samtools view -@ 12 -Sb ${output_pathway}/${sample_name}_reAligned.out.sam > ${output_pathway}/${sample_name}_re.bam
    rm ${output_pathway}/${sample_name}_reAligned.out.sam

    ### bedtools for re-mapped reads
    alignment_number_3=$(samtools view -c ${output_pathway}/${sample_name}_re.bam)
    bedtools bamtobed -split -cigar -bed12 -i ${output_pathway}/${sample_name}_re.bam > ${output_pathway}/${sample_name}_temp_re.bed
    awk '{print "chr" $0}' ${output_pathway}/${sample_name}_temp_re.bed > ${output_pathway}/${sample_name}_re.bed

    bedtools intersect -wa -wb -s -split \
        -a ${output_pathway}/${sample_name}_re.bed \
        -b /data/gladyshev/sz1012/damage_pipeline/rmsk_hg38.bed \
        > ${output_pathway}/${sample_name}_rmsk_temp.bed

    uniq ${output_pathway}/${sample_name}_rmsk_temp.bed > ${output_pathway}/${sample_name}_re_filtered.bed

    ### uniqmap processing
    alignment_number=$(samtools view -c ${output_pathway}/${sample_name}_uniqmap.bam)
    alignment_number_2=$(samtools view -c ${output_pathway}/${sample_name}Aligned.sortedByCoord.out.bam)

    bedtools bamtobed -split -cigar -bed12 \
        -i ${output_pathway}/${sample_name}_uniqmap.bam \
        > ${output_pathway}/${sample_name}_uniqmap_.bed

    awk '{print "chr" $0}' ${output_pathway}/${sample_name}_uniqmap_.bed > ${output_pathway}/${sample_name}_uniqmap.bed
    rm ${output_pathway}/${sample_name}_uniqmap_.bed

    bedtools intersect -wa -wb -s -split \
        -a ${output_pathway}/${sample_name}_uniqmap.bed \
        -b /data/gladyshev/sz1012/damage_pipeline/hg38_intron_stopsite_f.bed \
        > ${output_pathway}/${sample_name}_ss_temp.bed

    uniq ${output_pathway}/${sample_name}_ss_temp.bed > ${output_pathway}/${sample_name}_ss.bed

    python /data/gladyshev/sz1012/damage_pipeline/re2.py \
        -i ${output_pathway}/${sample_name}_re_filtered.bed \
        -o ${output_pathway}/${sample_name}_re.txt -l ${info_file}

    python /data/gladyshev/sz1012/damage_pipeline/ss2.py \
        -i ${output_pathway}/${sample_name}_ss.bed \
        -o ${output_pathway}/${sample_name}_ss.txt -l ${info_file}

    ### Arriba fusion detection
    /data/gladyshev/sz1012/damage_pipeline/arriba_v2.5.1/arriba \
        -x ${output_pathway}/${sample_name}Aligned.sortedByCoord.out.bam \
        -o ${output_pathway}/${sample_name}_fusions.tsv \
        -O ${output_pathway}/${sample_name}_fusions.discarded.tsv \
        -a /data/gladyshev/sz1012/damage_pipeline/arriba_v2.5.1/GRCh38.fa \
        -g /data/gladyshev/sz1012/damage_pipeline/arriba_v2.5.1/RefSeq_hg38.gtf \
        -b /data/gladyshev/sz1012/damage_pipeline/arriba_v2.5.1/database/blacklist_hg38_GRCh38_v2.5.1.tsv.gz \
        -k /data/gladyshev/sz1012/damage_pipeline/arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
        -t /data/gladyshev/sz1012/damage_pipeline/arriba_v2.5.1/database/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
        -p /data/gladyshev/sz1012/damage_pipeline/arriba_v2.5.1/database/protein_domains_hg38_GRCh38_v2.5.1.gff3

    ### RSEM expression
    rsem-calculate-expression --paired-end --star -p 25  ${file_pathway}/${sample_name}_1.fastq  ${file_pathway}/${sample_name}_2.fastq  ${star_reference_pathway2}  ${output_pathway}/${sample_name}

    ### Junction analysis
    python /data/gladyshev/sz1012/damage_pipeline/convert_sj_to_bed.py \
        ${output_pathway}/${sample_name}SJ.out.tab \
        ${output_pathway}/${sample_name}SJ.bed

    bedtools intersect -a ${output_pathway}/${sample_name}SJ.bed \
                       -b /data/gladyshev/sz1012/damage_pipeline/human_exons.bed -wa \
                       > ${output_pathway}/${sample_name}_temp1.bed

    uniq ${output_pathway}/${sample_name}_temp1.bed > ${output_pathway}/${sample_name}_temp2.bed

    awk '{print "chr" $0}' ${output_pathway}/${sample_name}_temp2.bed > ${output_pathway}/${sample_name}SJ_exon_uniq_.bed

    bedtools intersect -a ${output_pathway}/${sample_name}SJ_exon_uniq_.bed \
                       -b /data/gladyshev/sz1012/damage_pipeline/UP000005640_9606_domain_sort.bed \
                       -wa -wb > ${output_pathway}/${sample_name}_temp3.bed

    bedtools intersect -a ${output_pathway}/${sample_name}SJ_exon_uniq_.bed \
                       -b /data/gladyshev/sz1012/damage_pipeline/UP000005640_9606_domain_sort.bed \
                       -wa | uniq > ${output_pathway}/${sample_name}SJ_exon_domain.bed

    awk '!($2 < $8 && $3 > $9)' ${output_pathway}/${sample_name}_temp3.bed | sort -u \
        > ${output_pathway}/${sample_name}SJ_exon_domain_uniq.bed

    ### Summary
    n1=$(awk '{sum+=$5} END{print sum}' ${output_pathway}/${sample_name}SJ_exon_domain_uniq.bed)
    n2=$(awk '{sum+=$5} END{print sum}' ${output_pathway}/${sample_name}SJ_exon_domain.bed)
    n3=$(awk '{sum+=$5} END{print sum}' ${output_pathway}/${sample_name}SJ.bed)

    echo ${sample_name} number ${alignment_number} ${alignment_number_2} ${alignment_number_3} ${n1} ${n2} ${n3} >> ${info_file}

    ### CLEAN UP
    rm ${output_pathway}/${sample_name}_uniqmap.bam
    rm ${output_pathway}/${sample_name}_temp_re.bed
    rm ${output_pathway}/${sample_name}_rmsk_temp.bed
    rm ${output_pathway}/${sample_name}_ss_temp.bed
    rm ${output_pathway}/${sample_name}_temp1.bed
    rm ${output_pathway}/${sample_name}_temp2.bed
    rm ${output_pathway}/${sample_name}_temp3.bed
    rm ${output_pathway}/${sample_name}.transcript.bam
    rm ${output_pathway}/${sample_name}_re.bam

    echo "${sample_name} is done"
done


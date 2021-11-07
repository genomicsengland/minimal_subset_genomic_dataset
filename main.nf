#!/usr/bin/env nextflow

// Define Channels from input
Channel
    .fromPath(params.table_vcf_location)
    .ifEmpty { exit 1, "Cannot find input file : ${params.table_vcf_location}" }
    .splitCsv(skip:1)
    .map {file_name, vcf, vcf_idx -> [ file_name, file(vcf), file(vcf_idx) ] }
    .take( params.number_of_files_to_process )
    .set { ch_input_list }

Channel
    .fromPath(params.region_file_location)
    .ifEmpty { exit 1, "Cannot find input file : ${params.region_file_location}" }
    .set { ch_region_file }

// Define Process
process copy_vcf {
    tag "$file_name"
    publishDir "${params.outdir}/copy", mode: 'copy'

    input:
    set val(file_name), file(vcf), file(vcf_idx) from ch_input_list
    
    output:
    set val(file_name), file("${file_name}_copy.vcf.gz"), file("${file_name}_copy.vcf.gz.csi") into ch_copied
    
    script:
    """
    cp ${vcf} ${file_name}_copy.vcf.gz
    cp ${vcf_idx} ${file_name}_copy.vcf.gz.csi
    """
  }

process subset_vcf {
    tag "$file_name"
    publishDir "${params.outdir}/subset", mode: 'copy'

    input:
    set val(file_name), file(vcf), file(vcf_idx) from ch_copied
    each file(region_file) from ch_region_file

    output:
    set file("${file_name}_subset.vcf.gz"), file("${file_name}_subset.vcf.gz.csi") into ch_out

    script:
    """
    bcftools view -R ${region_file} ${vcf} -Oz -o ${file_name}_subset.vcf.gz
    bcftools index ${file_name}_subset.vcf.gz
    """
  }


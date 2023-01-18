#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// define workflow
workflow {
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
    
    if (params.skip_reformat) {
        SUBSET(ch_input_list, ch_region_file)
    }
    else {
        REFORMAT(ch_input_list)
        SUBSET(REFORMAT.out, ch_region_file)
    }
}

// Define Processes
process REFORMAT {
    tag "$file_name"
    publishDir "${params.outdir}/reformat", mode: 'copy'
        
    input:
    tuple val(file_name), file(vcf), file(vcf_idx)
        
    output:
    tuple val(file_name), file("${file_name}_reformat.bcf.gz"), file("${file_name}_reformat.bcf.gz.csi")
        
    script:
    """
    bcftools view ${vcf} -Ob -o ${file_name}_reformat.bcf.gz
    bcftools index ${file_name}_reformat.bcf.gz
    rm -f ${vcf} ${vcf_idx}
    """
}


process SUBSET {
    tag "$file_name"
    publishDir "${params.outdir}/subset", mode: 'copy'
    
    input:
    tuple val(file_name), file(input_file), file(input_file_idx)
    each file(region_file)
    
    output:
    tuple file("${file_name}_subset.vcf.gz"), file("${file_name}_subset.vcf.gz.csi")
    
    script:
    """
    bcftools view -T ${region_file} ${input_file} -Oz -o ${file_name}_subset.vcf.gz
    bcftools index ${file_name}_subset.vcf.gz
    rm -f ${input_file} ${input_file_idx} ${region_file}
    """
}

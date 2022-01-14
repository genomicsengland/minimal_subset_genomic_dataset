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

// Define Processes
if (params.skip_reformat) {
    ch_reformatted = ch_input_list
}
else {
    process reformat {
        tag "$file_name"
        publishDir "${params.outdir}/reformat", mode: 'copy'
        
        input:
        set val(file_name), file(vcf), file(vcf_idx) from ch_input_list
        
        output:
        set val(file_name), file("${file_name}_reformat.bcf.gz"), file("${file_name}_reformat.bcf.gz.csi") into ch_reformatted
        
        script:
        """
        bcftools view ${vcf} -Ob -o ${file_name}_reformat.bcf.gz
        bcftools index ${file_name}_reformat.bcf.gz
        rm -f ${vcf} ${vcf_idx}
        """
    }
}

process subset {
    tag "$file_name"
    publishDir "${params.outdir}/subset", mode: 'copy'
    
    input:
    set val(file_name), file(input_file), file(input_file_idx) from ch_reformatted
    each file(region_file) from ch_region_file
    
    output:
    set file("${file_name}_subset.vcf.gz"), file("${file_name}_subset.vcf.gz.csi") into ch_out
    
    script:
    """
    bcftools view -T ${region_file} ${input_file} -Oz -o ${file_name}_subset.vcf.gz
    bcftools index ${file_name}_subset.vcf.gz
    rm -f ${input_file} ${input_file_idx} ${region_file}
    """
}

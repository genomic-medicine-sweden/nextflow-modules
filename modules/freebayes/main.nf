//
// Initialize options with default values.
//
def initParams(Map params) {
    params.args = params.args ?: ''
    params.publishDir = params.publishDir ?: ''
    params.publishDirMode = params.publishDirMode ?: ''
    params.publishDirOverwrite = params.publishDirMode ?: false
    return params
}

params = initParams(params)

process freebayes {
  label "process_medium"
  tag "${sampleName}"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(fasta)
    tuple path(bam), path(bai) 

  output:
    tuple val(sampleName), path(output)

  script:
    output = "${sampleName}.vcf"
    """
    freebayes ${params.args.join(' ')} -f ${fasta} ${bam} > ${output}
    """
}

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
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(fasta)
    tuple path(bam), path(bai) 

  output:
    tuple val(sampleName), path(output)

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}.vcf"
    """
    freebayes ${args} -f ${fasta} ${bam} > ${output}
    """
}

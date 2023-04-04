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

process post_align_qc {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(bam)
    path bai
    path reference

  output:
    tuple val(sampleName), path(output)

  script:
    output = "${sampleName}_bwa.qc"
    """
    postaln_qc.pl ${bam} ${reference} ${sampleName} ${task.cpus} > ${output}
    """
}
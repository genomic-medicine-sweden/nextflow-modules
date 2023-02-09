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

process skesa {
  tag "${sampleName}"
  label "process_high"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)

  output:
    tuple val(sampleName), path("${output}")

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "${reads[0]},${reads[1]}" : "${reads[0]}"
    output = "${sampleName}.fa"
    """
    skesa --reads ${inputData} ${args} > ${output}
    """
}
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

process kraken {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    path database

  output:
    tuple val(sampleName), path("${output}"), emit: output
    tuple val(sampleName), path("${report}"), emit: report

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}_kraken.out"
    report = "${sampleName}_kraken.report"
    """
    kraken2 \\
    ${args} \\
    --threads ${task.cpus} \\
    --db ${database} \\
    --output ${output} \\
    --report ${report} \\
    ${reads.join(' ')}
    """
}

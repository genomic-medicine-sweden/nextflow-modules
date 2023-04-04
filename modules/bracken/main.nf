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

process bracken {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(report)
    path database

  output:
    tuple val(sampleName), path("${output}"), emit: output
    tuple val(sampleName), path("${outputReport}"), emit: report

  script:
    output = "${sampleName}_bracken.out"
    outputReport = "${sampleName}_bracken.report"
    """
    bracken \\
    ${args.join(' ')} \\
    -d ${database} \\
    -i ${report} \\
    -o ${output} \\
    -w ${outputReport}
    """
}

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
  label "process_high"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path report
    path database

  output:
    path("${output}"), emit: output
    path("${outputReport}"), emit: report

  script:
    sampleName = report.simpleName
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

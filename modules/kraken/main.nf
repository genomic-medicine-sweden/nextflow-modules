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
  label "process_high"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    path database

  output:
    path("${output}")
    path("${report}")

  script:
    output = "${sampleName}_kraken.out"
    report = "${sampleName}_kraken.report"
    """
    kraken2 \\
    ${args.join(' ')} \\
    --threads ${task.cpus} \\
    --db ${database} \\
    --output ${output} \\
    --report ${report} \\
    ${reads.join(' ')}
    """
}

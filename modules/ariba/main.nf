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

process ariba_prepareref {
  tag "${fasta.simpleName}"
  label "process_high"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path fasta
    path metadata

  output:
    path(outputDir)

  script:
    def outputDir = params.outdir ? params.outdir : "ariba_reference"
    def metadata = metadata ? "--metadata ${metadata}" : "--all_coding"
    """
    ariba prepareref \\
    --threads ${task.cpus} \\
    ${args.join(' ')} \\
    ${metadata} \\
    --fasta ${fasta} \\
    ${outputDir}
    """
}

process ariba_run {
  tag "${sampleName}"
  label "process_high"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    path referenceDir

  output:
    path("${sampleName}_ariba_report.tsv")

  script:
    outputName = params.outdir ? params.outdir : "ariba_output"
    """
    ariba run ${args.join(' ')} --threads ${task.cpus} ${referenceDir} ${reads.join(' ')} ${outputName}
    cp ${outputName}/report.tsv ${sampleName}_ariba_report.tsv
    """
}

process ariba_summary {
  tag "${report.simpleName}"
  label "process_high"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path report

  output:
    path("${outputPrefix}.csv")

  script:
    outputPrefix = params.prefix ?: report.simpleName.replaceFirst('report', 'summary')
    """
    ariba summary ${outputPrefix} ${report}
    """
}

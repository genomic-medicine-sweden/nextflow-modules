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
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path fasta
    path metadata

  output:
    path(outputDir)

  script:
    def outputDir = "ariba_reference"
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
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    path referenceDir

  output:
    tuple val(sampleName), path("${sampleName}_ariba_report.tsv")

  script:
    outputName = "ariba_output"
    """
    ariba run ${args.join(' ')} --force --threads ${task.cpus} --verbose ${referenceDir} ${reads.join(' ')} ${outputName}
    cp ${outputName}/report.tsv ${sampleName}_ariba_report.tsv
    """
}

process ariba_summary {
  tag "${sampleName}"
  label "process_high"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(report)

  output:
    tuple val(sampleName), path("${outputPrefix}.csv")

  script:
    outputPrefix = params.prefix ?: report.simpleName.replaceFirst('report', 'summary')
    """
    ariba summary ${args.join(' ')} ${outputPrefix} ${report}
    """
}

process ariba_summary_to_json {
  tag "${sampleName}"
  label "process_low"
  publishDir params.publishDir, 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(report), path(summary) 
    path reference

  output:
    tuple val(sampleName), path("${output}"), emit: output

  script:
    output = "${summary.simpleName}_export.json"
    """
    ariba2json.pl ${reference} ${summary} ${report} > ${output}
    """
}
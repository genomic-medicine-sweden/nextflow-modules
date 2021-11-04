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

process quast {
  tag "${assembly.simpleName}"
  label "process_medium"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path assembly
    path reference

  output:
    path "${assembly.simpleName}.quast.tsv"

  script:
    sampleName = "${assembly.simpleName}"
    args = params.args ? params.args : ''
    reference = reference ? "-r ${reference}" : ''
    outputDir = 'quast_outdir'
    """
    quast.py ${args.join(' ')} ${assembly} ${reference} -o ${outputDir} -t ${task.cpus}
    cp ${outputDir}/transposed_report.tsv ${sampleName}.quast.tsv
    """
}

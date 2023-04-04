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
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(assembly)
    path reference

  output:
    tuple val(sampleName), path("${assembly.simpleName}.quast.tsv")

  script:
    def args = task.ext.args ?: ''
    sampleName = "${assembly.simpleName}"
    reference = reference ? "-r ${reference}" : ''
    outputDir = 'quast_outdir'
    """
    quast.py ${args} ${assembly} ${reference} -o ${outputDir} -t ${task.cpus}
    cp ${outputDir}/transposed_report.tsv ${sampleName}.quast.tsv
    """
}

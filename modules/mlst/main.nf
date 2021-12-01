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

process mlst {
  tag "${sampleName}"
  label "process_medium"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(assembly)
    val specie

  output:
    tuple val(sampleName), path('*.tsv'), optional: true, emit: tsv
    tuple val(sampleName), path('*.json'), optional: true, emit: json
    tuple val(sampleName), path('*.novel'), optional: true, emit: novel

  script:
    outputName = "${assembly.simpleName}.mlst"
    blastDbPath = params.blastdb ? "--blastdb ${params.blastdb}" : ""
    pubmlstDataDir = params.pubmlstData ? "--datadir ${params.datadir}" : ""
    """
    mlst \\
    ${params.args.join(' ')} \\
    ${blastDbPath} \\
    ${pubmlstDataDir} \\
    --scheme ${specie} \\
    --json ${outputName}.json \\
    --novel ${outputName}.novel \\
    --threads ${task.cpus} \\
    ${assembly}
    """
}

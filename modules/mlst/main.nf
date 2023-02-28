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

def getAbbrevSpeciesName(fullName) {
  "Convert the full name to the abbreviated version"
  names = fullName.split(' ')
  return names[0][0] + names[1]
}

params = initParams(params)

process mlst {
  tag "${sampleName}"
  label "process_medium"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(assembly)
    val species
    path blastDb

  output:
    tuple val(sampleName), path('*.tsv'), optional: true, emit: tsv
    tuple val(sampleName), path('*.json'), optional: true, emit: json
    tuple val(sampleName), path('*.novel'), optional: true, emit: novel

  script:
    outputName = "${sampleName}.mlst"
    abbrevName = getAbbrevSpeciesName(species)
    blastDbPath = blastDb ? "--blastdb ${blastDb}/mlst.fa" : ""
    //pubmlstDataDir = pubmlstData ? "--datadir ${pubmlstData}" : ""
    //${pubmlstDataDir} \\
    """
    mlst \\
      ${params.args.join(' ')} \\
      ${blastDbPath} \\
      --scheme  ${abbrevName} \\
      --json ${outputName}.json \\
      --novel ${outputName}.novel \\
      --threads ${task.cpus} \\
      ${assembly}
    """
}

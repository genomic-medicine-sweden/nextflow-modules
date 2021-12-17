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

process virulencefinder {
  label "process_medium"
  tag "${sampleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    val databases
    path virulenceDb

  output:
    tuple val(sampleName), path(outputFile)
    
  script:
    def databasesArgs = databases ? "--databases ${databases.join(',')}" : ""
    def outputDir = "results"
    outputFile = "virulencefinder_${sampleName}.json"
    """
    mkdir ${outputDir}
    virulencefinder.py              \\
    --infile ${reads.join(' ')}     \\
    ${databasesArgs}                \\
    --databasePath ${virulenceDb}   \\
    --outputPath ${outputDir}
    cp ${outputDir}/data.json ${outputFile}
    """
}

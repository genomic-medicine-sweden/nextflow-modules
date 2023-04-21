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
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    val databases
    path virulenceDb

  output:
    tuple val(sampleName), path(outputFile), emit: json
    tuple val(sampleName), path(metaFile)  , emit: meta
    
  script:
    databasesArgs = databases ? "--databases ${databases.join(',')}" : ""
    outputFile = "virulencefinder_${sampleName}.json"
    metaFile = "virulencefinder_meta_${sampleName}.json"
    """
    # Get db version
    DB_HASH=\$(git -C ${virulenceDb} rev-parse HEAD)
    JSON_FMT='{"name": "%s", "version": "%s", "type": "%s"}'
    printf "\$JSON_FMT" "virulencefinder" "\$DB_HASH" "database" > ${metaFile}

    # Run virulencefinder
    virulencefinder.py              \\
    --infile ${reads.join(' ')}     \\
    ${databasesArgs}                \\
    --databasePath ${virulenceDb}
    cp data.json ${outputFile}
    """

 stub:
    outputFile = "virulencefinder_${sampleName}.json"
    metaFile = "virulencefinder_meta_${sampleName}.json"
    """
    touch $outputFile
    touch $metaFile
    """
}

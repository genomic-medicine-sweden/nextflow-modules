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

process resfinder {
  label "process_medium"
  tag "${sampleName}"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    val specie
    path resfinderDb
    path pointfinderDb

  output:
    tuple val(sampleName), path(outputFileJson), emit: json
    tuple val(sampleName), path(metaFile), emit: meta
    path outputFileGene, emit: geneTable
    path outputFilePoint, emit: pointTable
    
  script:
    def resfinderFinderParams = pointfinderDb ? "--acquired --db_path_res ${resfinderDb}" : ""
    def pointFinderParams = pointfinderDb ? "--point --db_path_point ${pointfinderDb}" : ""
    def specieArgs = specie ? "--species '${specie}'" : ""
    outputFileJson = "resfinder_${sampleName}.json"
    metaFile = "resfinder_meta_${sampleName}.json"
    outputFileGene = "pheno_table_${sampleName}.txt"
    outputFilePoint = "point_table_${sampleName}.txt"
    """
    # Get db version
    RES_DB_HASH=\$(git -C ${resfinderDb} rev-parse HEAD)
    POINT_DB_HASH=\$(git -C ${pointfinderDb} rev-parse HEAD)
    JSON_FMT='[{"name": "%s", "version": "%s", "type": "%s"},{"name": "%s", "version": "%s", "type": "%s"}]'
    printf "\$JSON_FMT" "resfinder" \$RES_DB_HASH "database" "pointfinder" \$POINT_DB_HASH "database" > $metaFile
    # Get resfinder path
    RESF=\$(which run_resfinder.py)

    # Run resfinder
    \$RESF                          \\
    --inputfastq ${reads.join(' ')} \\
    ${specieArgs}                   \\
    ${resfinderFinderParams}        \\
    ${pointFinderParams}
    cp std_format_under_development.json ${outputFileJson}
    cp pheno_table.txt ${outputFileGene}
    cp PointFinder_results.txt ${outputFilePoint}
    """
}

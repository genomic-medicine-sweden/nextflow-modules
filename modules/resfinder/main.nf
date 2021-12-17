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
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    val specie
    path resfinderDb
    path pointfinderDb

  output:
    tuple val(sampleName), path("${outputPath}/std_format_under_development.json"), emit: json
    path "${outputPath}/pheno_table.txt", emit: geneTable
    path "${outputPath}/PointFinder_results.txt", emit: pointTable
    
  script:
    def pointFinderParams = pointfinderDb ? "--point --db_path_point ${pointfinderDb}" : ""
    def specieArgs = specie ? "--species ${specie}" : ""
    outputPath = "results"
    """
    run_resfinder.py                \\
    --inputfastq ${reads.join(' ')} \\
    ${specieArgs}                   \\
    --db_path_res ${resfinderDb}    \\
    ${pointFinderParams}            \\
    --outputPath ${outputPath}
    """
}

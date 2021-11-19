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

process chewbbaca_allelecall {
  label "process_medium"
  tag "${assembly.simpleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path input
    path schemaDir
    path trainingFile

  output:
    path('output_dir/*/results_alleles.tsv'), emit: results
    path("${missingLoci}"), emit: missing

  script:
    id = "${assembly.simpleName}"
    output = "${id}.vcf"
    missingLoci = "chewbbaca.missingloci"
    trainingFile = $trainingFile ?: "--ptf ${trainingFile}"
    """
    echo ${input} > batch_input.list
    flock -e ${params.localTempDir}/chewbbaca.lock \\
      chewBBACA.py AlleleCall \\
      ${params.join(' ')} \\
      --cpu ${task.cpus} \\
      --output-directory output_dir \\
      ${trainingFile} \\
      --schema-directory ${schemaDir}
    bash parse_missing_loci.sh batch_input.list 'output_dir/*/results_alleles.tsv' ${missingloci}
    """
}

process chewbbaca_split_results {
  label "process_low"
  tag "${assembly.simpleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path input

  output:
    path("${output}")

  script:
    id = "${input.simpleName}"
    output = "${id}.chewbbaca"
    """
    head -1 ${input} > ${output}
    grep ${id} ${input} >> ${output}
    """
}

process chewbbaca_split_missing_loci {
  label "process_low"
  tag "${assembly.simpleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path input

  output:
    path("${output}")

  script:
    id = "${input.simpleName}"
    output = "${id}.chewbbaca"
    """
    grep ${id} ${input} > ${output}
    """
}

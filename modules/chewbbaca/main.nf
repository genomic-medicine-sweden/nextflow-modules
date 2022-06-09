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
  tag "${sampleName}"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(input)
    path schemaDir
    path trainingFile

  output:
    tuple val(sampleName), path('output_dir/*/results_alleles.tsv')
    //tuple val(sampleName), path("${missingLoci}"), emit: missing

  script:
    output = "${sampleName}.vcf"
    missingLoci = "chewbbaca.missingloci"
    trainingFile = trainingFile ? "--ptf ${trainingFile}" : ""
    """
    echo ${input} > batch_input.list
    flock -e ${params.localTempDir}/chewbbaca.lock \\
      chewBBACA.py AlleleCall \\
      -i batch_input.list \\
      ${params.args.join(' ')} \\
      --cpu ${task.cpus} \\
      --output-directory output_dir \\
      ${trainingFile} \\
      --schema-directory ${schemaDir}
    #bash parse_missing_loci.sh batch_input.list 'output_dir/*/results_alleles.tsv' ${missingLoci}
    """
}

process chewbbaca_split_results {
  label "process_low"
  tag "${assembly.simpleName}"
  publishDir "${params.publishDir}", 
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
  publishDir "${params.publishDir}", 
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

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

process bwa_index {
  label "process_high"
  tag "${reference.simpleName}"

  input:
    path reference

  output:
    path "${reference}.*", emit: reference

  script:
    """
    bwa index ${reference}
    """
}

process bwa_mem {
  label "process_high"
  tag "${sampleName}"

  input:
    tuple val(sampleName), path(reads)
    path referenceIdx

  output:
    path "${sampleName}.sam"

  script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
      ${args.join(' ')} \\
      -t ${task.cpus} \\
      \${INDEX} \\
      ${reads.join(' ')} > ${sampleName}.sam 
    """
}

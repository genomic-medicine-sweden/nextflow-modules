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
  tag "${sampleName}"

  input:
    tuple val(sampleName), path(reference)

  output:
    tuple val(sampleName), path("${reference}.*"), emit: reference

  script:
    """
    bwa index ${reference} ${reference.baseName}/${reference}
    """
}

process bwa_mem {
  label "process_high"
  tag "${sampleName}"

  input:
    tuple val(sampleName), path(reads)
    path referenceIdx

  output:
    tuple val(sampleName), path("${sampleName}.sam")

  script:
    def args = task.ext.args ?: ''
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
      ${args} \\
      -t ${task.cpus} \\
      \${INDEX} \\
      ${reads.join(' ')} > ${sampleName}.sam 
    """
}

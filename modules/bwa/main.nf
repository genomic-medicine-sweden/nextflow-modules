process index {
  label "process_high"
  tag "${reference.simpleName}"

  input:
    path reference

  output:
    path "${reference}.*"

  script:
    """
    bwa index ${reference}
    """
}

process mem {
  label "process_high"
  tag "${sampleName}"

  input:
    path referenceDir
    tuple val(sampleName), path(reads)

  output:
    path "${samleName}.sam"

  script:
    // define names and paths
    referenceFile = ''
    readOne = reads[0]
    readTwo = reads[1]
    // compute flags
    markShort = params.markShort ? '-M' : ''
    """
    bwa mem ${markShort} -t ${task.cpus} ${referenceFile} ${readOne} ${readTwo} > ${sampleName}.sam
    """

}

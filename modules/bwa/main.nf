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

process BWAMEM {
  label "process_high"
  tag "${sampleName}"

  input:
    tuple val(sampleName), path(r1), path(r2)

  output:
    path "${sampleName}.sam"

  script:
    // define names and paths
    //referenceFile = ''
    //readOne = reads[0]
    //readTwo = reads[1]
    // compute flags
    markShort = params.markShort ? '-M' : ''

    """
    bwa mem ${markShort} -t ${task.cpus} ${params.genome} ${r1} ${r2} > ${sampleName}.sam 
    """

}

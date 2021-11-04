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

process bwaMem {
  label "process_high"
  tag "${sampleName}"

  input:
    tuple val(sampleName), path(r1), path(r2)

  output:
    path "${sampleName}.sam"

  script:

    markShort = params.markShort ? '-M' : ''

    """
    bwa mem ${markShort} -t ${task.cpus} ${params.genome} ${r1} ${r2} > ${sampleName}.sam 
    """

}

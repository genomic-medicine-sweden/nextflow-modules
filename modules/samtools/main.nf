// Samtools functions

process view {
  tag "$input"
  label "process_medium"

  input:
    path input
    path fasta

  output:
    path('*.bam'), optional: true, emit: bam
    path('*.cram'), optional: true, emit: cram

  script:
    def reference = fasta ? "--reference ${fasta} -C" : ""
    def prefix = input.simpleName
    def fileType = input.getExtension()
    """
    samtools view $reference ${input} > ${prefix}.${fileType}
    """
}

process sort {
  tag "$input"
  label "process_medium"

  input:
    path input
    path fasta

  output:
    path('*.bam'), optional: true, emit: bam
    path('*.cram'), optional: true, emit: cram

  script:
    def reference = fasta ? "--reference ${fasta} -C" : ""
    def prefix = input.simpleName
    def fileType = fasta ? "cram" : input.getExtension()
    def maxMemory = params.maxMemory ? "-m ${params.maxMemory}" : ""
    """
    samtools sort ${reference} -@ $task.cpus -o ${prefix}.sorted.${fileType} ${input}
    """
}

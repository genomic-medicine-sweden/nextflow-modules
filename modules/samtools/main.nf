// Samtools functions

process samtools_view {
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

process samtools_sort {
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
    """
    samtools sort ${reference} -@ $task.cpus -o ${prefix}.sorted.${fileType} ${input}
    """
}

process samtools_index {
  tag "$input"
  label "process_medium"

  input:
    path input

  output:
    path output

  script:
    output = "${input}.bai"
    """
    samtools index -@ $task.cpus ${input}
    """
}

process samtools_faidx {
  tag "$input"
  label "process_low"

  input:
    path input

  output:
    path output

  script:
    output = "${input}.fai"
    """
    samtools faidx ${input}
    """
}

// Samtools functions

process samtools_view {
  tag "$input"
  scratch params.scratch

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
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(input)
    path fasta

  output:
    tuple val(sampleName), path('*.bam'), optional: true, emit: bam
    tuple val(sampleName), path('*.cram'), optional: true, emit: cram

  script:
    def reference = fasta ? "--reference ${fasta} -O cram" : "-O bam"
    def prefix = input.simpleName
    def fileType = fasta ? "cram" : "bam"
    """
    samtools sort ${reference} -@ $task.cpus -o ${prefix}.sorted.${fileType} ${input}
    """
}

process samtools_index {
  tag "${sampleName}"
  scratch params.scratch

  input:
    tuple val(sampleName), path(input)

  output:
    tuple val(sampleName), path(output)

  script:
    output = "${input}.bai"
    """
    samtools index -@ $task.cpus ${input}
    """
}

process samtools_faidx {
  tag "$input"
  scratch params.scratch

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

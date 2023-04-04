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

process spades_iontorrent {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads), val(platform)

  output:
    tuple val(sampleName), path("${sampleName}.fasta")

  when:
    platform == "iontorrent"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    outputDir = params.publishDir ? params.publishDir : 'spades'
    """
    spades.py ${args} ${inputData} -t ${task.cpus} -o ${outputDir}
    mv ${outputDir}/contigs.fasta ${sampleName}.fasta
    """
}

process spades_illumina {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads), val(platform)

  output:
    tuple val(sampleName), path("${sampleName}.fasta")

  when:
    task.ext.when && platform == "illumina"

  script:
    def args = task.ext.args ?: ''
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    outputDir = params.publishDir ? params.publishDir : 'spades'
    """
    spades.py ${args} ${inputData} -t ${task.cpus} -o ${outputDir}
    mv ${outputDir}/contigs.fasta ${sampleName}.fasta
    """
}


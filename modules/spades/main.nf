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

process spades {
  tag "${sampleName}"
  label "process_high"
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads)
    val platform

  output:
    tuple val(sampleName), path("${sampleName}.fasta")

  script:
    def args = task.ext.args ?: ''
    def platform_command = platform ? "--${platform}" : ""
    def inputData = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    outputDir = params.publishDir ? params.publishDir : 'spades'
    """
    spades.py ${platform_command} ${args.join(' ')} ${inputData} -t ${task.cpus} -o ${outputDir}
    mv ${outputDir}/contigs.fasta ${sampleName}.fasta
    """
}

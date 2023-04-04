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

process assembly_trim_clean {
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads), val(platform)

  output:
    tuple val(sampleName), path(output), val(platform)

  when:
    platform == "iontorrent"

  script:
    def args = task.ext.args ?: ''
    output = "${sampleName}_cleaned.fastq.gz"
    """
    run_assembly_trimClean.pl --numcpus ${task.cpus} ${args} -i ${reads} -o ${output}
    """
}
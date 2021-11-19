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

process freebayes {
  label "process_medium"
  tag "${assembly.simpleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path assembly
    path assemblyIdx
    path mappedReads

  output:
    path output

  script:
    id = "${assembly.simpleName}"
    output = "${id}.vcf"
    """
    freebayes ${params.join(' ')} -f ${assembly} ${mappedReads} > ${output}
    """
}

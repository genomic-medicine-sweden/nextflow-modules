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


def getFileBasename(path) {
  "Return the basename of a filepath"
  return path.split('/')[-1] - ~/\.\w+$/
}

params = initParams(params)

process mask_polymorph_assembly {
  label "process_low"
  tag "${assembly.simpleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path assembly
    path polymorph

  output:
    path output

  script:
    id = "${assembly.simpleName}"
    output = "${id}_masked.fasta"
    """
    error_corr_assembly.pl ${assembly} ${polymorph} > ${output}
    """
}

process export_to_cdm {
  label "process_low"
  tag "${assembly.simpleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path cgmlstMissingLoci
    path quast
    path postQc

  output:
    path output

  script:
    id = "${assembly.simpleName}"
    output = "${id}.cdm"
    rundir = 'fool'

    """
    echo --run-folder ${rundir} \\
         --sample-id ${id} \\
         --assay microbiology \\
         --qc ${postQc} \\
         --asmqc ${quast} \\
         --micmisloc ${cgmlstMissingLoci} > ${output}
    """
}

process export_to_cgviz {
  label "process_low"
  tag "${quast}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    //path cgmlst
    path quast
    path mlst
    path bracken
    path cgmlst
    path ariba

  output:
    path output

  script:
    id = "${quast}"
    output = "${id}.cgviz"
    rundir = 'fool'
    """
    echo --overwrite \\
         --sample-id ${id} \\
         --species ${params.specie} \\
         --run ${rundir} \\
         --in ${cgmlst} \\
         --kraken ${bracken} \\
         --aribavir ${ariba} \\
         --micmisloc ${cgmlstMissingLoci} \\
         --quast ${quast} > ${output}
    """ 
}

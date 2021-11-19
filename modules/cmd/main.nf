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

process ariba_summary_to_json {
  tag "${summary.simpleName}"
  label "process_low"
  publishDir "${params.outdir}",
    mode: params.publishDirMode,
    overwrite: params.publishDirOverwrite

  input:
    path summary
    path report
    path reference

  output:
    path("${output}"), emit: output

  script:
    output = "${summary.simpleName}_export.json"
    """
    ariba2json.pl ${reference} ${summary} ${report} > ${output}
    """
}

process post_align_qc {
  tag "${bam.simpleName}"
  label "process_low"
  publishDir "${params.outdir}",
    mode: params.publishDirMode,
    overwrite: params.publishDirOverwrite

  input:
    path bam
    path index
    path cgmlst

  output:
    path report

  script:
    id = "${bam.simpleName}"
    output = "${id}_bwa.qc"
    """
    postaln_qc.pl ${bam} ${reference} ${id} ${task.cpus} > ${output}
    """
}

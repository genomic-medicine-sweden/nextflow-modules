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
  tag "${sampleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(assembly), path(polymorph)

  output:
    tuple val(sampleName), path(output)

  script:
    output = "${sampleName}_masked.fasta"
    """
    error_corr_assembly.pl ${assembly} ${polymorph} > ${output}
    """
}

process export_to_cdm {
  label "process_low"
  tag "${sampleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(cgmlstMissingLoci), path(quast), path(postQc)

  output:
    path(output)

  script:
    output = "${sampleName}.cdm"
    rundir = 'fool'

    """
    echo --run-folder ${rundir} \\
         --sample-id ${sampleName} \\
         --assay microbiology \\
         --qc ${postQc} \\
         --asmqc ${quast} \\
         --micmisloc ${cgmlstMissingLoci} > ${output}
    """
}

process export_to_cgviz {
  label "process_low"
  tag "${sampleName}"
  publishDir "${params.outdir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    //path cgmlst
    tuple val(sampleName), path(quast), path(mlst), path(cgmlst), path(virulence), path(resistance)
    //path bracken

  output:
    path(output)

  script:
    output = "${sampleName}.cgviz"
    rundir = 'fool'
    //--kraken ${bracken} \\
    //--micmisloc ${cgmlstMissingLoci} \\
    """
    echo --overwrite \\
         --sample-id ${sampleName} \\
         --species ${params.specie} \\
         --run ${rundir} \\
         --cgmlst ${cgmlst} \\
         --virulence ${virulence} \\
         --resistance ${resfinder} \\
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

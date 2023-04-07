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
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
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
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
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
  tag "${sampleName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    path runInfo
    file meta
    //paths
    tuple val(sampleName), val(quast), val(mlst), val(cgmlst), val(virulence), val(resistance)

  output:
    path(output)

  script:
    output = "${sampleName}_cgviz.json"
    //--kraken ${bracken} \\
    quastArgs = quast ? "--quast ${quast}" : "" 
    mlstArgs = mlst ? "--mlst ${mlst}" : "" 
    cgmlstArgs = cgmlst ? "--cgmlst ${cgmlst}" : "" 
    resfinderArgs = resistance ? "--resistance ${resistance}" : "" 
    virulenceArgs = virulence ? "--virulence ${virulence}" : "" 
    metaArgs = meta ? "--process-metadata  ${meta[1..-1].join(' --process-metadata ')}" : ""
    """
    combine_results.py \\
      --sample-id ${sampleName} \\
      --run-metadata ${runInfo} \\
      ${metaArgs} \\
      ${quastArgs} \\
      ${mlstArgs} \\
      ${cgmlstArgs} \\
      ${virulenceArgs} \\
      ${resfinderArgs} \\
      ${output}
    """ 
}

process ariba_summary_to_json {
  tag "${summary.simpleName}"
  scratch params.scratch
  publishDir "${params.publishDir}",
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
  scratch params.scratch
  publishDir "${params.publishDir}",
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

process save_analysis_metadata {
  tag "${workflow.runName}"
  scratch params.scratch
  publishDir "${params.publishDir}", 
    mode: params.publishDirMode, 
    overwrite: params.publishDirOverwrite

  input:
    tuple val(sampleName), path(reads), val(platform)  

  output:
    path(output)

  script:
    def seqType = reads.size() == 2 ? "PE" : "SE"
    output = "${sampleName}_analysis_meta.json"
    """
    #!/usr/bin/env python
    import json

    res = {
        "workflow_name": "$workflow.runName",
        "sample_name": "${sampleName}",
        "sequencing_platform": "${platform}",
        "sequencing_type": "${seqType}",
        "date": "$workflow.start",
        "pipeline": "$workflow.scriptName",
        "version": "$workflow.manifest.version",
        "commit": "$workflow.commitId",
        "configuration_files": "$workflow.configFiles"[1:-1].split(','),
        "analysis_profile": "$workflow.profile",
        "command": "$workflow.commandLine",
    }
    with open("$output", 'w') as out:
        json.dump(res, out, indent=2)
    """
}

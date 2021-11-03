/*
Testing scripts for modules created with the test data in two servers - trannel (testing) and hopper (server). See nextflow.config associated with this script to find the configuration of different servers and settings
*/

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */ 

params.genome = '/fs1/resources/ref/micro/species/saureus/ref.fasta'
//params.markShort = '-M'



/* 
 * Import modules 
 */
include { BWAMEM } from '../../modules/bwa/main.nf'

workflow {
  Channel
        .fromPath(params.csv)
        .splitCsv(header:true)
        .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
        .set { input }

      // Align with the reference genome 

      BWAMEM(input)
}

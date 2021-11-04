// Example of an bacterial analysis pipeline
nextflow.enable.dsl=2

// define workflow parameters
params.outdir = 'test_output'
params.specie = 'ecoli'

include { samtools_sort } from './modules/samtools/main.nf'
include { spades } from './modules/spades/main.nf' addParams( args: ['--iontorrent', '--only-assembler'] )
include { quast } from './modules/quast/main.nf' addParams( args: [] )
include { mlst } from './modules/mlst/main.nf' addParams( args: [] )
include { ariba_run } from './modules/ariba/main.nf' addParams( outdir: 'ariba_test_outdir', args: ['--force'] )
include { ariba_summary } from './modules/ariba/main.nf' addParams( args: ['--col_filter', 'n', '--row_filter', 'n'] )
include { kraken } from './modules/kraken/main.nf' addParams( args: ['--gzip-compressed'] )


workflow bacterial_pipeline {
  //bamFile = channel.fromPath('MT21000025.dedup.bam', checkIfExists: true)
  read_pairs_ch = Channel .fromFilePairs('SRR*_{1,2}.fastq.gz').view()
  assembly = channel.fromPath('SRR10490537.fasta', checkIfExists: true)
  genomeReference = channel.fromPath('/fs1/resources/ref/micro/species/saureus/ref.fasta', checkIfExists: true)
  aribaReferenceDir = channel.fromPath('/fs1/resources/ref/micro/species/saureus/ariba', checkIfExists: true)
  krakenDb = channel.fromPath('/fs1/resources/ref/micro/krakenstd', checkIfExists: true)


  main:
    //samtools_sort(bamFile, [])
    //assembly = spades(read_pairs_ch)
    //quast(assembly, reference)
    //mlst(assembly, params.specie)
    //ariba_summary(ariba_run(read_pairs_ch, aribaReferenceDir))
    kraken(read_pairs_ch, krakenDb)
}

## Pipelines at CMD 

These folder consists of the test data and pipelines used to test the modules created. The basic and elementary nextflow run inside a pipeline is as such :

```
  nextflow run main.nf  --csv samplesheet.csv -profile standard --genome /fs1/resources/ref/micro/species/saureus/ref.fasta
```
# nextflow-modules
Shared nextflow modules and assets used at CMD

# The basic structure of the pipeline is as follows
```

├── LICENSE
├── modules
│   ├── bwa
│   │   └── main.nf
│   ├── samtools
│   │   └── main.nf
│   └── senteion
│       └── bwa
├── pipeline
│   ├── micro
│   │   ├── data
│   │   │   └── micro
│   │   │       ├── SRR10490537_1.fastq.gz
│   │   │       └── SRR10490537_2.fastq.gz
│   │   ├── main.nf
│   │   └── nextflow.config
│   └── nextflow.config
└── README.md

```
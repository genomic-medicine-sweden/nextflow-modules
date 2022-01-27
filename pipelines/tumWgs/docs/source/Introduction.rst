Introduction
====================================================================
This pipeline is developed using nextflow modular system and consists of different subworkflows. The main workflow consists mainly of the sentieon distributed alignment system followed by tnscope and dnnascoppe variant calling. following the bam/cram files emitted from the pipelines, we then call SNV variant subworkflow . Furtjer this worklow conists of annotation and filteraation with panel of normals.
Sunbsequenlty CNV and fustion variant are called upon.

Here is the general outline that has been generalized from `Reproducible research`_ at NBIs/ELIXIR 

.. image:: reproducibility-overview.png
   :width: 600


How to use the pipeline ?
====================================================================

Data
----------------------------------------------------------------------
The data are the fasta files and are in the seqdata files.n The parent directory for the seqdata for example at CMD is at the following site normally:: 
   
   /fs1/seqdata/


Environment
----------------------------------------------------------------------
All the softwares requried for the provided in the singualrity container. The defination fies for creating a singulaity container are present in the container folder. The container files at CMD are at following folder::

   /fs1/resources/container/

If noen of the container is present, build from the scratch using the container foledr in the 

Codes
----------------------------------------------------------------------
Source codes for running the pipelines and all necessary scripts are present in the files. The main.nf files consist of worflow for running the pipelines. The script for running pipelines and maintainaing scripts are present at folder::
   
    /fs1/bjorn/bnf-scripts/

Results
----------------------------------------------------------------------
The results from the pipelines are presented in the CoYOte . The simlink for the bam files, vcf files are presented in the result folder. The results from different analysis and pipelines are presenet in the folloeing folder::
   
   /fs1/results/


.. _Reproducible research: https://nbis-reproducible-research.readthedocs.io/en/course_2104/
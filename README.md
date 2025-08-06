## Microbiome Transcriptomic Integration 

This workflow pipeline does the inegration of metataxonomic and transcriptomics data 

## Installation 

1. Install Nextflow 
2. Run main as nexflow run main.nf 
3. Wait for execution. 
=======
1. Install conda environment microtrans through yaml microtrans.yml
2. Activate the microtrans environment 
3. Run modules using nextflow run *name of module*
4. Wait for execution.

## Run 

# Make bin dir executable
chmod +x bin/*
# Pre-create conda envirorments
mamba env create -f envs/microtrans.yml --name microtrans
# Activate microtrans envirorment 
mamba activate microtrans
# Run
nextflow run workflows/main.nf -c config/nextflow.config 


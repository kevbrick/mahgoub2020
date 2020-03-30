# Mahgoub et al. 2020 (zcwpw1) nextflow pipeline : 
This nextflow pipeline is configured to work on a SLURM-based HPC with modules. It can be relatively easily configured to run on other systems (see nextflow documentation : https://www.nextflow.io/). 

## Requirements:
*Tools / programs :*
- bedtools	2.27.1
- deeptools	3.0.1
- macs	2.1.2
- meme	5.0.1
- minimap2	2.17
- nextflow	0.30.2
- picard	2.9.2
- R	3.5.2
- samtools	1.9
- sratoolkit	2.9.2
- ucsc	388

**NOTE:** These are the recommended (and tested) version of each tool, however newer / older versions will work in most cases. One exception is older versions of nextflow, which will not work. 

## Global variables required: 
$SLURM_JOBID   : Specifies the temporary subfolder to use  (see Temp folder requirements below)

### Temp folder requirements
The pipeline requires a high-level temporary folder called /lscratch. On a SLURM-based HPC, each job is assigned a global id ($SLURM_JOBID) and this is appended to the temp folder name for each process. This is currently hard-coded. Thus, there is a requirement for :

/lscratch folder for temporary files
SLURM_JOBID global variable for each HPC job.

## Pipeline execution

### RUN ON SLURM CLUSTER
```
nextflow run -c accessoryFiles/config/nextflow.config \
    Zcwpw1_CutNRunPipe.nf \
    --outdir output \
    --getData \
    --prdm9BG accessoryFiles/B6_PRDM9ChIPSeq.bedgraph \
    --projectdir `pwd` \
    -with-report Zcwpw1Pipeline.nxfReport.html \
    -with-trace Zcwpw1Pipeline.nxfTrace.html \
    -with-timeline Zcwpw1Pipeline.nxfTimeline.html
```

### RUN ON LOCAL MACHINE
```
nextflow run -c accessoryFiles/config/nextflow.local.config \
    Zcwpw1_CutNRunPipe.nf \
    --outdir output \
     --getData \
    --prdm9BG accessoryFiles/B6_PRDM9ChIPSeq.bedgraph \
    --projectdir `pwd` \
    -with-report Zcwpw1Pipeline.nxfReport.html \
    -with-trace Zcwpw1Pipeline.nxfTrace.html \
    -with-timeline Zcwpw1Pipeline.nxfTimeline.html
```







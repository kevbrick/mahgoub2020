#!/bin/bash

module load nextflow

nextflow run -c accessoryFiles/nxfConfig/nextflow.config Zcwpw1_CutNRunPipe.nf  --outdir output --getData  --projectdir `pwd`  -with-report Zcwpw1Pipeline.nxfReport.html  -with-trace Zcwpw1Pipeline.nxfTrace.html -with-timeline Zcwpw1Pipeline.nxfTimeline.html

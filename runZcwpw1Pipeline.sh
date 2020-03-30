#!/bin/bash
module load nextflow/0.30.2
nextflow run -resume -c accessoryFiles/nxfConfig/nextflow.config Zcwpw1_CutNRunPipe.nf --outdir output --prdm9BG accessoryFiles/B6_PRDM9ChIPSeq.bedgraph --projectdir `pwd` -with-report Zcwpw1Pipeline.nxfReport.html -with-trace Zcwpw1Pipeline.nxfTrace.html -with-timeline Zcwpw1Pipeline.nxfTimeline.html -with-dag Zcwpw1_pipe.dot

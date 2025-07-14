#!/bin/bash
SINGULARITY=`which singularity`
sudo ionice -c 3 $SINGULARITY build -F yahs_pipeline.sif Singularity

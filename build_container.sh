#!/bin/bash
SINGULARITY=`which singularity`
sudo $SINGULARITY build -F yahs_pipeline.sif Singularity

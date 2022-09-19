#!/bin/bash

SCRIPTDIR="/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/"
LOGDIR="/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/log"
SCRIPT="run_counting.sh"

qsub \
-N cellrangerarcsamp123457 \
-M carlos.ramirez@bioquant.uni-heidelberg.de \
-m bea \
-l walltime=72:00:00,mem=600g \
-l nodes=1:ppn=32 \
-o $LOGDIR \
-e $LOGDIR \
$SCRIPTDIR$SCRIPT 

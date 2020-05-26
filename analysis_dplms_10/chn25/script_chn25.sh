#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd_dplms/processPSD /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_10/chn25 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms_10.tier.root 25 1 10 > /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_10/chn25/psd_chn25.out

chmod 666 *.*


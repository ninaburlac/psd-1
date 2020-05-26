#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd_dplms/processPSD /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_100/chn16 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms_100.tier.root 16 1 100 > /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_100/chn16/psd_chn16.out

chmod 666 *.*


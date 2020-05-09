#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/dandrea/Analysis/PSD/processPSD /nfs/gerda6/users/dandrea/Analysis/PSD/analysis_dplms_10/chn8 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms_10.tier.root 8 1 10 > /nfs/gerda6/users/dandrea/Analysis/PSD/analysis_dplms_10/chn8/psd_chn8.out

chmod 666 *.*


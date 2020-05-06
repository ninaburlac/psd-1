#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/dandrea/Analysis/PSD/processPSD /nfs/gerda6/users/dandrea/Analysis/PSD/analysis_dplms/chn2 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms.tier.root 2 1 > /nfs/gerda6/users/dandrea/Analysis/PSD/analysis_dplms/chn2/psd_chn2.out

chmod 666 *.*


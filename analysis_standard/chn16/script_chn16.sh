#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/dandrea/Analysis/PSD/processPSD /nfs/gerda6/users/dandrea/Analysis/PSD/analysis_standard/chn16 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms.tier.root 16 0 > /nfs/gerda6/users/dandrea/Analysis/PSD/analysis_standard/chn16/psd_chn16.out

chmod 666 *.*


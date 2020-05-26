#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd_dplms/processPSD /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_sigma_1000/chn2 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms_1000.tier.root 2 1 1000 > /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_sigma_1000/chn2/psd_chn2.out

chmod 666 *.*


#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd_dplms/processPSD /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_sigma_1000/chn27 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms_1000.tier.root 27 1 1000 > /nfs/gerda6/users/burlac/psd_dplms/analysis_dplms_sigma_1000/chn27/psd_chn27.out

chmod 666 *.*


#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd_dplms/processPSD /nfs/gerda6/users/burlac/psd_dplms/analysis_standard_sigma/chn37 /nfs/gerda6/users/burlac/psd/analysis/run0098-cal_dplms_10.tier.root 37 0 0 > /nfs/gerda6/users/burlac/psd_dplms/analysis_standard_sigma/chn37/psd_chn37.out

chmod 666 *.*


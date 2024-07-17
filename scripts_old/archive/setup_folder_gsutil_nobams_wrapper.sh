release=$1
study_id=$2

workingfolder=~/my_lustre/analysis/${release}_${study_id}
gsutil_path=gs://vo_afun_release/${release}
logfolder=$workingfolder/folder_setup_log 
mkdir -p $logfolder 
~/my_lustre/cnv_calls/CNV_funestus/scripts/setup_folder_gsutil_nobams.sh $workingfolder $gsutil_path $study_id > $logfolder/setup_folder_gsutil.log 2>&1 & 

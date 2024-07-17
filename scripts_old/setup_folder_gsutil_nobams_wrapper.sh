release=$1
study_id=$2

workingfolder=~/personal/${release}_${study_id}
gsutil_path=gs://vo_afun_release/${release}
logfolder=$workingfolder/folder_setup_log 
mkdir -p $logfolder 
~/scripts/CNV_funestus/scripts/setup_folder_gsutil_nobams.sh $workingfolder $gsutil_path $study_id > $logfolder/setup_folder_gsutil.log 2>&1 & 

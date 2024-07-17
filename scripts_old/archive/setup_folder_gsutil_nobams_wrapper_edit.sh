≈release=$1
study_id=$2
sampleset=$3

workingfolder=~/my_lustre/analysis/${release}_${study_id}
gsutil_path=gs://vo_afun_release/${release}
logfolder=$workingfolder/folder_setup_log 
mkdir -p $logfolder 
~/my_lustre/cnv_calls/CNV_funestus/scripts/archive/setup_folder_gsutil_nobams_edit.sh $workingfolder $gsutil_path $study_id $sampleset > $logfolder/setup_folder_gsutil.log 2>&1 & 

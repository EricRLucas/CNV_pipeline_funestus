workingfolder=$1
manifest=$2
accessibility_file=$3
mapqprop_file=$4
GC_content_file=$5
sample_group_id=$6
scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts

source activate cnv37 

# Now calculate median coverage by GC 
echo "Calculating median and variance of coverage by GC bin for sample group ${sample_group_id}"
python $scriptsfolder/calculate_median_coverage_by_GC.py 0.9 \
                                                         $accessibility_file \
                                                         0.5 \
                                                         $mapqprop_file \
                                                         $manifest \
                                                         $GC_content_file \
                                                         $workingfolder \
                                                         $sample_group_id \
                                                         > $workingfolder/calculate_mean_coverage_by_GC_09_05_${sample_group_id}.log 2>&1



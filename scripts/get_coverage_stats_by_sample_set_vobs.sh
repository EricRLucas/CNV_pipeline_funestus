workingfolder=$1
manifest=$2
accessibility_file=$3
mapqprop_file=$4
GC_content_file=$5
scriptsfolder=~/scripts/CNV_funestus/scripts

source activate cnv37 

# Now calculate median coverage by GC 
echo "Calculating median and variance of coverage by GC bin"
python $scriptsfolder/calculate_median_coverage_by_GC.py 0.9 \
                                                         $accessibility_file \
                                                         0.5 \
                                                         $mapqprop_file \
                                                         $manifest \
                                                         $GC_content_file \
                                                         $workingfolder \
                                                         all \
                                                         > $workingfolder/calculate_mean_coverage_by_GC_09_05.log 2>&1



coveragefolder=$1
manifest=$2
coverage_variance_file=$3
metadata=$4
outputfolder=$5
gene_coordinates_file=$6

scriptsfolder=~/scripts/CNV_funestus/scripts

R-3.6.1 --version

R-3.6.1 --slave -f $scriptsfolder/modal_cnv.r --args $manifest \
                                                     $gene_coordinates_file \
                                                     $metadata \
                                                     $output_name \
                                                     $coverage_variance_file \
                                                     $coveragefolder \
                                                     $outputfolder \
                                                     > $outputfolder/modal_CNV.log 2>&1


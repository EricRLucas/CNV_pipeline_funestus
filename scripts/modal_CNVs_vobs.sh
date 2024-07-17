coveragefolder=$1
manifest=$2
coverage_variance_file=$3
metadata=$4
outputfolder=$5
gene_coordinates_file=$6

scriptsfolder=~/scripts/CNV_funestus/scripts

module load /software/spack_environments/default/00/share/spack/modules/linux-ubuntu22.04-x86_64_v3/r-3.6.3/python-3.11.6

R --version

R --slave -f $scriptsfolder/modal_cnv.r --args $manifest \
                                               $gene_coordinates_file \
                                               $metadata \
                                               $output_name \
                                               $coverage_variance_file \
                                               $coveragefolder \
                                               $outputfolder \
                                               > $outputfolder/modal_CNV.log 2>&1


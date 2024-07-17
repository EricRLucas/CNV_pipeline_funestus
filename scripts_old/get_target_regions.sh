workingfolder=$1
manifest=$2
coverage_variance_file=$3
gene_coordinates_file=$4
metadata=$5
ncores=$6

coveragefolder=$workingfolder/coverage
diagnostic_reads_folder=$workingfolder/diagnostic_reads
scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts
plotting_functions_file=$scriptsfolder/plotting_functions.r

cd $workingfolder

mkdir -p target_regions_analysis

export R_LIBS_USER="/nfs/users/nfs_e/el10/R-modules:$R_LIBS_USER"

R-3.6.1 --version

R-3.6.1 --slave -f $scriptsfolder/get_target_regions.r --args $manifest \
                                                              $gene_coordinates_file \
                                                              $metadata \
                                                              $coverage_variance_file \
                                                              $coveragefolder \
                                                              $diagnostic_reads_folder \
                                                              $plotting_functions_file \
                                                              $ncores \
                                                              > target_regions_analysis/get_target_regions.log 2>&1



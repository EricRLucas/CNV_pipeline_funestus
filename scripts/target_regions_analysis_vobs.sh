workingfolder=$1
manifest=$2
coverage_variance_file=$3
gene_coordinates_file=$4
metadata=$5
ncores=$6

coveragefolder=$workingfolder/coverage
diagnostic_reads_folder=$workingfolder/diagnostic_reads
scriptsfolder=~/scripts/CNV_funestus/scripts
plotting_functions_file=$scriptsfolder/plotting_functions.r

cd $workingfolder

mkdir -p target_regions_analysis

export R_LIBS_USER="~/R-modules:$R_LIBS_USER"

module load /software/spack_environments/default/00/share/spack/modules/linux-ubuntu22.04-x86_64_v3/r-3.6.3/python-3.11.6

R --version

R --slave -f $scriptsfolder/target_regions_analysis.r --args $manifest \
                                                             $gene_coordinates_file \
                                                             $metadata \
                                                             $coverage_variance_file \
                                                             $coveragefolder \
                                                             $diagnostic_reads_folder \
                                                             $plotting_functions_file \
                                                             $ncores \
                                                             > target_regions_analysis/target_regions_analysis.log 2>&1



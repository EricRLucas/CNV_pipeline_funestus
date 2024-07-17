coveragefolder=$1
manifest=$2
chrom=$3
output_name=$4
coverage_variance_file=$5
ncores=$6
metadata=$7
gene_coordinates_file=$8
detox_genes_file=$9

workingfolder=$coveragefolder/$chrom/HMM_output
outputfolder=$workingfolder/$output_name
scriptsfolder=~/scripts/CNV_funestus/scripts

export R_LIBS_USER="~/R-modules:$R_LIBS_USER"

module load /software/spack_environments/default/00/share/spack/modules/linux-ubuntu22.04-x86_64_v3/r-3.6.3/python-3.11.6

R --version

R --slave -f $scriptsfolder/CNV_analysis.r --args $chrom \
                                                  $manifest \
                                                  $coverage_variance_file \
                                                  $gene_coordinates_file \
                                                  $detox_genes_file \
                                                  $workingfolder \
                                                  $ncores \
                                                  $outputfolder \
                                                  $metadata \
                                                  > $coveragefolder/$chrom/CNV_analysis_logs/CNV_analysis_${output_name}.log 2>&1


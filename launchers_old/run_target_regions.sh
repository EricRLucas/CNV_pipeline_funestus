samplelist=$1
metadata_file=$2
study=$3

scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts
rootfolder=/nfs/users/nfs_j/jn11/my_lustre/analysis
workingfolder=$rootfolder/$study
diagnostic_reads_folder=$workingfolder/diagnostic_reads
coverage_variance_file=$workingfolder/coverage/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=$rootfolder/genome_data/Afun_gene_regions.csv
ncores=20
logfolder=$workingfolder/logfolders/get_target_regions
errorfolder=$workingfolder/errorfolders/get_target_regions

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/target_regions_output_%J.txt \
     -e $errorfolder/target_regions_error_%J.txt \
     -q normal \
     -n $ncores \
     -R"select[mem>30000] rusage[mem=30000] span[hosts=1]" \
     -M30000 \
     $scriptsfolder/get_target_regions.sh $workingfolder \
                                          $samplelist \
                                          $coverage_variance_file \
                                          $gene_coordinates_file \
                                          $metadata_file \
                                          $ncores

samplelist=$1
metadata_file=$2
study=$3

scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts
rootfolder=/nfs/users/nfs_j/jn11/my_lustre/analysis
workingfolder=$rootfolder/$study
diagnostic_reads_folder=$workingfolder/diagnostic_reads
coverage_variance_file=$workingfolder/coverage/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=$rootfolder/genome_data/Afun_gene_regions.csv
ncores=5
logfolder=$diagnostic_reads_folder/logfolders/target_regions_analysis
errorfolder=$diagnostic_reads_folder/errorfolders/target_regions_analysis

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/target_regions_analysis_output_%J.txt \
     -e $errorfolder/target_regions_analysis_error_%J.txt \
     -q long \
     -n $ncores \
     -R"select[mem>10000] rusage[mem=10000] span[hosts=1]" \
     -M10000 \
     $scriptsfolder/target_regions_analysis_vobs.sh $workingfolder \
                                                    $samplelist \
                                                    $coverage_variance_file \
                                                    $gene_coordinates_file \
                                                    $metadata_file \
                                                    $ncores

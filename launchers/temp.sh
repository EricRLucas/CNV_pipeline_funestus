samplelist=$1
metadata_file=$2
study=$3

scriptsfolder=~/scripts/CNV_funestus/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10/funestus
workingfolder=$rootfolder/$study
diagnostic_reads_folder=$workingfolder/diagnostic_reads
coverage_variance_file=$workingfolder/coverage/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=/lustre/scratch126/tol/teams/lawniczak/users/am60/vobs/cnv_production/analysis/genome_data/Afun_gene_regions.csv
ncores=12
logfolder=$diagnostic_reads_folder/logfolders/target_regions_analysis
errorfolder=$diagnostic_reads_folder/errorfolders/target_regions_analysis

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/target_regions_analysis_output_%J.txt \
     -e $errorfolder/target_regions_analysis_error_%J.txt \
     -q long \
     -n $ncores \
     -R"select[mem>50000] rusage[mem=50000] span[hosts=1]" \
     -M50000 \
     $scriptsfolder/target_regions_analysis_vobs.sh $workingfolder \
                                                    $samplelist \
                                                    $coverage_variance_file \
                                                    $gene_coordinates_file \
                                                    $metadata_file \
                                                    $ncores

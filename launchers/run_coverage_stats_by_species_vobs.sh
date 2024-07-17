study=$1

rootfolder=/lustre/scratch126/gsu/team112/personal/el10/funestus
logfolder=$rootfolder/$study/coverage/logfolders/coverage_stats
errorfolder=$rootfolder/$study/coverage/errorfolders/coverage_stats
workingfolder=$rootfolder/$study/coverage
manifestfolder=$rootfolder/$study/data
datafolder=/lustre/scratch126/tol/teams/lawniczak/users/am60/vobs/cnv_production/analysis/genome_data
scriptsfolder=~/scripts/CNV_funestus/scripts

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/coverage_stats_log.txt \
	 -e $errorfolder/coverage_stat_error.txt \
	 -R"select[mem>500] rusage[mem=500] span[hosts=1]" \
	 -M500 \
	 ${scriptsfolder}/get_coverage_stats_by_sample_set_vobs.sh $workingfolder \
                                                               $manifestfolder/sample_manifest.txt \
                                                               $datafolder/mean_accessibility.csv \
                                                               $datafolder/mapq_proportions_allchrom.csv \
                                                               $datafolder/Afun_genome_GC_content.csv \
                                                               all


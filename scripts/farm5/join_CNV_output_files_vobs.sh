study=$1

rootfolder=/lustre/scratch126/gsu/team112/personal/el10/funestus
coveragefolder=$rootfolder/$study/coverage
outputfolder=$rootfolder/$study/coverage_CNVs

mkdir -p $outputfolder

combined_full_CNV_table=$outputfolder/full_coverage_CNV_table.csv
combined_full_raw_CNV_table=$outputfolder/full_raw_CNV_table.csv

cat $coveragefolder/2RL/HMM_output/funestus_CNV/full_coverage_CNV_table_2RL.csv > $combined_full_CNV_table
tail -n +2 $coveragefolder/3RL/HMM_output/funestus_CNV/full_coverage_CNV_table_3RL.csv >> $combined_full_CNV_table
tail -n +2 $coveragefolder/X/HMM_output/funestus_CNV/full_coverage_CNV_table_X.csv >> $combined_full_CNV_table

cat $coveragefolder/2RL/HMM_output/funestus_CNV/full_raw_CNV_table_2RL.csv > $combined_full_raw_CNV_table
tail -n +2 $coveragefolder/3RL/HMM_output/funestus_CNV/full_raw_CNV_table_3RL.csv >> $combined_full_raw_CNV_table
tail -n +2 $coveragefolder/X/HMM_output/funestus_CNV/full_raw_CNV_table_X.csv >> $combined_full_raw_CNV_table

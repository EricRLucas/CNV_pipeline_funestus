#path to download bamfiles
bamfilepaths=$1
#number of jobs 
joblimit=$2
#study name: starts  with v sth
study=$3

#pipe tpo scripts
scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts
#pipe to output folder
rootfolder=/nfs/users/nfs_j/jn11/my_lustre/analysis
outputfolder=$rootfolder/$study
bamfilefolder=$outputfolder/bamlinks
coveragefolder=$outputfolder/coverage
SSFAfolder=$outputfolder/diagnostic_reads/SSFA
breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
logfolder=$outputfolder/logfolders/coverage_and_diagnostic_reads
errorfolder=$outputfolder/errorfolders/coverage_and_diagnostic_reads
allchrom=(2RL 3RL X)

echo "Calculating coverage and identifying diagnostic reads for samples from ${bamfilepaths} on `date`." >> $outputfolder/added_samples.log

# Create the output folders if necessary
mkdir -p $bamfilefolder
mkdir -p $logfolder
mkdir -p $errorfolder
for chrom in ${allchrom[@]}
do
	mkdir -p $coveragefolder/$chrom/coveragelogs
done
mkdir -p $SSFAfolder/2RL/Cyp6_region/SSFAlogs #looks at all discordant reads same strand, face away, cross chromosome
mkdir -p $breakpointsfolder/2RL/Cyp6_region/breakpointlogs  #soft clip reads
# Get the number of bamfiles that need processing

numbams=($(wc -l $bamfilepaths))
echo "This sample set contains ${numbams} bam files." >> $outputfolder/added_samples.log

bsub -J "bamArray[1-$numbams]%$joblimit" \
     -R"select[mem>500] rusage[mem=500]" \
     -M500 \
     -o $logfolder/log_%J.%I.txt \
     -e $errorfolder/error_%J.%I.txt \
     ' '${scriptsfolder}'/get_windowed_coverage_and_diagnostic_reads.sh '${bamfilefolder}' '${bamfilepaths}' ${LSB_JOBINDEX} '$outputfolder


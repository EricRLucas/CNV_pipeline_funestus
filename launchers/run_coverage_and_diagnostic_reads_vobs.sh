bamfilepaths=$1
joblimit=$2
study=$3

scriptsfolder=~/scripts/CNV_funestus/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10/funestus
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
mkdir -p $SSFAfolder/2RL/Cyp6_region/SSFAlogs
mkdir -p $SSFAfolder/2RL/Ache_region/SSFAlogs
mkdir -p $SSFAfolder/2RL/Esterase_region/SSFAlogs
mkdir -p $SSFAfolder/2RL/Gst1_region/SSFAlogs
mkdir -p $SSFAfolder/2RL/Harb1_region/SSFAlogs
mkdir -p $SSFAfolder/2RL/Gst2_region/SSFAlogs
mkdir -p $SSFAfolder/2RL/Zinccarbo_region/SSFAlogs
mkdir -p $SSFAfolder/3RL/GABA_region/SSFAlogs
mkdir -p $SSFAfolder/3RL/Carboxypep_region/SSFAlogs
mkdir -p $SSFAfolder/X/Cyp9_region/SSFAlogs
mkdir -p $SSFAfolder/X/RDGA_region/SSFAlogs
mkdir -p $SSFAfolder/X/NADHCyp_region/SSFAlogs

mkdir -p $breakpointsfolder/2RL/Cyp6_region/breakpointlogs
mkdir -p $breakpointsfolder/2RL/Ache_region/breakpointlogs
mkdir -p $breakpointsfolder/2RL/Esterase_region/breakpointlogs
mkdir -p $breakpointsfolder/2RL/Gst1_region/breakpointlogs
mkdir -p $breakpointsfolder/2RL/Harb1_region/breakpointlogs
mkdir -p $breakpointsfolder/2RL/Gst2_region/breakpointlogs
mkdir -p $breakpointsfolder/2RL/Zinccarbo_region/breakpointlogs
mkdir -p $breakpointsfolder/3RL/GABA_region/breakpointlogs
mkdir -p $breakpointsfolder/3RL/Carboxypep_region/breakpointlogs
mkdir -p $breakpointsfolder/X/Cyp9_region/breakpointlogs
mkdir -p $breakpointsfolder/X/RDGA_region/breakpointlogs
mkdir -p $breakpointsfolder/X/NADHCyp_region/breakpointlogs

# Get the number of bamfiles that need processing
numbams=($(wc -l $bamfilepaths))
echo "This sample set contains ${numbams} bam files." >> $outputfolder/added_samples.log

bsub -J "bamArray[1-$numbams]%$joblimit" \
     -R"select[mem>10000] rusage[mem=10000]" \
     -M10000 \
     -o $logfolder/log_%J.%I.txt \
     -e $errorfolder/error_%J.%I.txt \
     ' '${scriptsfolder}'/get_windowed_coverage_and_diagnostic_reads.sh '${bamfilefolder}' '${bamfilepaths}' ${LSB_JOBINDEX} '$outputfolder


#submits jobs
#links jobs to sample names
 
bamfilefolder=$1
bamfilepaths=$2
samplenum=$3 
outputfolder=$4
scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts

samplename=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f1))
bampath=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f2))
bamfile=${bamfilefolder}/${samplename}.fixmate.bam

echo Downloading bam file for $samplename from $bampath to $bamfile

# Download the required bam: we dont need this
#curl ${bampath} > ${bamfile} 2> ${bamfilefolder}/${samplename}_download.log
#curl ${bampath}.bai > ${bamfile}.bai 2>> ${bamfilefolder}/${samplename}_download.log

#we need to create a conda environment for the next parts
source activate cnv37 

# Get the coverage data
echo Calculating coverage
allchrom=(2RL 3RL X)
coveragefolder=$outputfolder/coverage
for chrom in ${allchrom[@]}
do
	python ${scriptsfolder}/counts_for_HMM.py \
	       $bamfile \
	       $chrom \
	       300 300 10 \
	       ${coveragefolder}/${chrom}/counts_for_HMM_${samplename}_${chrom}_output.csv \
	       > ${coveragefolder}/${chrom}/coveragelogs/counts_for_HMM_${samplename}_${chrom}.log 2>&1
done

echo Identifying discordant reads
# Get the discordant reads
SSFA_script=$scriptsfolder/SSFA.py
SSFAfolder=$outputfolder/diagnostic_reads/SSFA
python $SSFA_script $bamfile 2RL 8300000:9200000 ${SSFAfolder}/2RL/Cyp6_region/${samplename}_Cyp6_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Cyp6_region/SSFAlogs/${samplename}_Cyp6_SSFA_output.log 2>&1

# Get the soft clipped reads
breakpoints_script=$scriptsfolder/breakpoint_detector.py
breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
python $breakpoints_script $bamfile 2RL 8300000:9200000 ${breakpointsfolder}/2RL/Cyp6_region/${samplename}_Cyp6_breakpoints_output 10 > ${breakpointsfolder}/2RL/Cyp6_region/breakpointlogs/${samplename}_Cyp6_breakpoints_output.log 2>&1

# Delete the bam file 
#echo Deleting bam file
#rm ${bamfile}*

#echo


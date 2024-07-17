bamfilefolder=$1
bamfilepaths=$2
samplenum=$3 
outputfolder=$4
scriptsfolder=~/scripts/CNV_funestus/scripts

samplename=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f1))
bampath=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f2))
bamfile=${bamfilefolder}/${samplename}.fixmate.bam

echo Downloading bam file for $samplename from $bampath to $bamfile

# Download the required bam
curl ${bampath} > ${bamfile} 2> ${bamfilefolder}/${samplename}_download.log
curl ${bampath}.bai > ${bamfile}.bai 2>> ${bamfilefolder}/${samplename}_download.log

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
python $SSFA_script $bamfile 2RL 8240000:9820000 ${SSFAfolder}/2RL/Cyp6_region/${samplename}_Cyp6_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Cyp6_region/SSFAlogs/${samplename}_Cyp6_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2RL 19400000:19500000 ${SSFAfolder}/2RL/Ache_region/${samplename}_Ache_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Ache_region/SSFAlogs/${samplename}_Ache_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2RL 39020000:39850000 ${SSFAfolder}/2RL/Esterase_region/${samplename}_Esterase_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Esterase_region/SSFAlogs/${samplename}_Esterase_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2RL 42000000:42460000 ${SSFAfolder}/2RL/Gst1_region/${samplename}_Gst1_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Gst1_region/SSFAlogs/${samplename}_Gst1_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2RL 57054400:58000000 ${SSFAfolder}/2RL/Harb1_region/${samplename}_Harb1_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Harb1_region/SSFAlogs/${samplename}_Harb1_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2RL 76000000:76900000 ${SSFAfolder}/2RL/Gst2_region/${samplename}_Gst2_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Gst2_region/SSFAlogs/${samplename}_Gst2_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2RL 85450000:85500000 ${SSFAfolder}/2RL/Zinccarbo_region/${samplename}_Zinccarbo_SSFA_output.csv 10 > ${SSFAfolder}/2RL/Zinccarbo_region/SSFAlogs/${samplename}_Zinccarbo_SSFA_output.log 2>&1
python $SSFA_script $bamfile 3RL 13400000:14200000 ${SSFAfolder}/3RL/GABA_region/${samplename}_GABA_SSFA_output.csv 10 > ${SSFAfolder}/3RL/GABA_region/SSFAlogs/${samplename}_GABA_SSFA_output.log 2>&1
python $SSFA_script $bamfile 3RL 30300000:30350000 ${SSFAfolder}/3RL/Carboxypep_region/${samplename}_Carboxypep_SSFA_output.csv 10 > ${SSFAfolder}/3RL/Carboxypep_region/SSFAlogs/${samplename}_Carboxypep_SSFA_output.log 2>&1
python $SSFA_script $bamfile X 8350000:8870000 ${SSFAfolder}/X/Cyp9_region/${samplename}_Cyp9_SSFA_output.csv 10 > ${SSFAfolder}/X/Cyp9_region/SSFAlogs/${samplename}_Cyp9_SSFA_output.log 2>&1
python $SSFA_script $bamfile X 13590000:14000000 ${SSFAfolder}/X/RDGA_region/${samplename}_RDGA_SSFA_output.csv 10 > ${SSFAfolder}/X/RDGA_region/SSFAlogs/${samplename}_RDGA_SSFA_output.log 2>&1
python $SSFA_script $bamfile X 14350000:14600000 ${SSFAfolder}/X/NADHCyp_region/${samplename}_NADHCyp_SSFA_output.csv 10 > ${SSFAfolder}/X/NADHCyp_region/SSFAlogs/${samplename}_NADHCyp_SSFA_output.log 2>&1

# Get the soft clipped reads
breakpoints_script=$scriptsfolder/breakpoint_detector.py
breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
python $breakpoints_script $bamfile 2RL 8240000:9820000 ${breakpointsfolder}/2RL/Cyp6_region/${samplename}_Cyp6_breakpoints_output 10 > ${breakpointsfolder}/2RL/Cyp6_region/breakpointlogs/${samplename}_Cyp6_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2RL 19400000:19500000 ${breakpointsfolder}/2RL/Ache_region/${samplename}_Ache_breakpoints_output 10 > ${breakpointsfolder}/2RL/Ache_region/breakpointlogs/${samplename}_Ache_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2RL 39020000:39850000 ${breakpointsfolder}/2RL/Esterase_region/${samplename}_Esterase_breakpoints_output 10 > ${breakpointsfolder}/2RL/Esterase_region/breakpointlogs/${samplename}_Esterase_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2RL 42000000:42460000 ${breakpointsfolder}/2RL/Gst1_region/${samplename}_Gst1_breakpoints_output 10 > ${breakpointsfolder}/2RL/Gst1_region/breakpointlogs/${samplename}_Gst1_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2RL 57054400:58000000 ${breakpointsfolder}/2RL/Harb1_region/${samplename}_Harb1_breakpoints_output 10 > ${breakpointsfolder}/2RL/Harb1_region/breakpointlogs/${samplename}_Harb1_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2RL 76000000:76900000 ${breakpointsfolder}/2RL/Gst2_region/${samplename}_Gst2_breakpoints_output 10 > ${breakpointsfolder}/2RL/Gst2_region/breakpointlogs/${samplename}_Gst2_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2RL 85450000:85500000 ${breakpointsfolder}/2RL/Zinccarbo_region/${samplename}_Zinccarbo_breakpoints_output 10 > ${breakpointsfolder}/2RL/Zinccarbo_region/breakpointlogs/${samplename}_Zinccarbo_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 3RL 13400000:14200000 ${breakpointsfolder}/3RL/GABA_region/${samplename}_GABA_breakpoints_output 10 > ${breakpointsfolder}/3RL/GABA_region/breakpointlogs/${samplename}_GABA_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 3RL 30300000:30350000 ${breakpointsfolder}/3RL/Carboxypep_region/${samplename}_Carboxypep_breakpoints_output 10 > ${breakpointsfolder}/3RL/Carboxypep_region/breakpointlogs/${samplename}_Carboxypep_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile X 8350000:8870000 ${breakpointsfolder}/X/Cyp9_region/${samplename}_Cyp9_breakpoints_output 10 > ${breakpointsfolder}/X/Cyp9_region/breakpointlogs/${samplename}_Cyp9_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile X 13590000:14000000 ${breakpointsfolder}/X/RDGA_region/${samplename}_RDGA_breakpoints_output 10 > ${breakpointsfolder}/X/RDGA_region/breakpointlogs/${samplename}_RDGA_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile X 14350000:14600000 ${breakpointsfolder}/X/NADHCyp_region/${samplename}_NADHCyp_breakpoints_output 10 > ${breakpointsfolder}/X/NADHCyp_region/breakpointlogs/${samplename}_NADHCyp_breakpoints_output.log 2>&1

# Delete the bam file 
echo Deleting bam file
rm ${bamfile}*



study=$1

scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts
rootfolder=/nfs/users/nfs_j/jn11/my_lustre/analysis
coveragefolder=$rootfolder/$study/coverage
manifestfolder=$rootfolder/$study/data
logfolder=$coveragefolder/logfolders/HMM
errorfolder=$coveragefolder/errorfolders/HMM

Afun_GC_file=$rootfolder/genome_data/Afun_genome_GC_content.csv

allchrom=(2RL 3RL X)

# Create the outputfolders if necessary
mkdir -p $logfolder
mkdir -p $errorfolder

sample_manifest=$manifestfolder/sample_manifest.txt
echo "Running HMM for samples from $sample_manifest on `date`." >> $coveragefolder/HMM_added_samples.log

# Get some of the necessary input filenames
coverage_by_GC_file=median_coverage_by_GC_masked_09_05_all.csv
coverage_variance_file=coverage_variance_masked_09_05_all.csv
mapq_prop_file=$rootfolder/genome_data/mapq_proportions_allchrom.csv

# Run the HMM on each chromosome
for chrom in ${allchrom[@]}
do
	mkdir -p $coveragefolder/$chrom/HMM_output
	mkdir -p $coveragefolder/$chrom/HMM_logs

	bsub -o $logfolder/HMM_output_%J.txt \
		 -e $errorfolder/HMM_error_%J.txt \
		 -R"select[mem>500] rusage[mem=500] span[hosts=1]" \
		 -M500 \
		 -q basement \
		 ${scriptsfolder}/coverage_HMM_vobs.sh $coveragefolder \
							   $sample_manifest \
							   $chrom \
							   $Afun_GC_file \
							   $coverage_by_GC_file \
							   $coverage_variance_file \
							   $mapq_prop_file
done


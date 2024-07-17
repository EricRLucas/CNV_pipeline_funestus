metadata_file=$1
study=$2

scriptsfolder=~/scripts/CNV_funestus/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10/funestus
coveragefolder=$rootfolder/$study/coverage
manifestfolder=$rootfolder/$study/data
coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv
datafolder=/lustre/scratch126/tol/teams/lawniczak/users/am60/vobs/cnv_production/analysis/genome_data
gene_coordinates_file=$datafolder/Afun_gene_regions.csv
detox_genes_file=$datafolder/Afun_detox_genes.txt
ncores=2
logfolder=$coveragefolder/logfolders/CNV_analysis
errorfolder=$coveragefolder/errorfolders/CNV_analysis

allchrom=(2RL 3RL X)

mkdir -p $logfolder
mkdir -p $errorfolder

sample_manifest=$manifestfolder/sample_manifest.txt

output_id=funestus_CNV

for chrom in ${allchrom[@]}
do
	mkdir -p $coveragefolder/$chrom/CNV_analysis_logs

	bsub -o $logfolder/CNV_analysis_output_${chrom}_%J.txt \
		 -e $errorfolder/CNV_analysis_error_${chrom}_%J.txt \
		 -n $ncores \
		 -q long \
		 -R"select[mem>500] rusage[mem=500] span[hosts=1]" \
		 -M500 \
		 ${scriptsfolder}/coverage_CNVs_vobs.sh $coveragefolder \
		                                        $sample_manifest \
                                                $chrom \
                                                $output_id \
							                    $coverage_variance_file \
							                    $ncores \
							                    $metadata_file \
							                    $gene_coordinates_file \
							                    $detox_genes_file
done

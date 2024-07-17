metadata_file=$1
study=$2

scriptsfolder=~/scripts/CNV_funestus/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10/funestus
coveragefolder=$rootfolder/$study/coverage
sample_manifest=$rootfolder/$study/data/sample_manifest.txt
coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=/lustre/scratch126/tol/teams/lawniczak/users/am60/vobs/cnv_production/analysis/genome_data/Afun_gene_regions.csv
outputfolder=$rootfolder/$study/modal_CNVs
logfolder=$coveragefolder/logfolders/modal_CNV
errorfolder=$coveragefolder/errorfolders/modal_CNV

mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $outputfolder

bsub -o $logfolder/modal_CNVs_%J.txt \
     -e $errorfolder/modal_CNVs_%J.txt \
     -q long \
     -R"select[mem>10000] rusage[mem=10000] span[hosts=1]" \
     -M10000 \
     ${scriptsfolder}/modal_CNVs_vobs.sh $coveragefolder \
	                                     $sample_manifest \
	                                     $coverage_variance_file \
	                                     $metadata_file \
	                                     $outputfolder \
	                                     $gene_coordinates_file


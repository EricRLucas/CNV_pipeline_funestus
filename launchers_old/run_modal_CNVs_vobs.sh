metadata_file=$1
study=$2

scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts
rootfolder=/nfs/users/nfs_j/jn11/my_lustre/analysis
coveragefolder=$rootfolder/$study/coverage
sample_manifest=$rootfolder/$study/data/sample_manifest.txt
coverage_variance_file=$coveragefolder/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=$rootfolder/genome_data/Afun_gene_regions.csv
outputfolder=$rootfolder/$study/modal_CNVs
logfolder=$coveragefolder/logfolders/modal_CNV
errorfolder=$coveragefolder/errorfolders/modal_CNV

mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $outputfolder

bsub -o $logfolder/modal_CNVs_%J.txt \
     -e $errorfolder/modal_CNVs_%J.txt \
     -q basement \
     -R"select[mem>10000] rusage[mem=10000] span[hosts=1]" \
     -M10000 \
     ${scriptsfolder}/modal_CNVs_vobs.sh $coveragefolder \
	                                     $sample_manifest \
	                                     $coverage_variance_file \
	                                     $metadata_file \
	                                     $outputfolder \
	                                     $gene_coordinates_file


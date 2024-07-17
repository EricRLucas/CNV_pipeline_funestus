coveragefolder=$1
manifest=$2
chrom=$3
GC_content_file=$4
coverage_by_GC_file=$coveragefolder/$5
coverage_variance_file=$coveragefolder/$6
mapq_prop_file=$7

scriptsfolder=/nfs/users/nfs_j/jn11/my_lustre/cnv_calls/CNV_funestus/scripts

source activate cnv37 

python ${scriptsfolder}/HMM_process.py \
       $manifest \
       $chrom \
       $coveragefolder \
       $GC_content_file \
       $coverage_by_GC_file \
       $coverage_variance_file \
       $mapq_prop_file \
       0.5 \
       > ${coveragefolder}/${chrom}/HMM_logs/HMM_${chrom}.log 2>&1


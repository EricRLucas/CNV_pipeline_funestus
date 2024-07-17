# This script prepares a folder in which a given sample set will be run
folder_location=$1
gsutil_path=$2
sampleset=$3

# turn on extended globbing
shopt -s extglob

# Create the folder
mkdir -p $folder_location/data
# Download the metadata to the folder
sample_metadata=$gsutil_path/metadata/general/$sampleset/samples.meta.csv
gsutil cp $sample_metadata $folder_location/data
metadata_filename="$(basename -- $sample_metadata)"
# Create a sample manifest from the metadata
sed -e '1d' -e '2,$s/,.*//' $folder_location/data/$metadata_filename > $folder_location/data/sample_manifest.txt

# Download the table linking sample names with bam files
sample_bampaths=$gsutil_path/metadata/general/$sampleset/wgs_snp_data.csv
gsutil cp $sample_bampaths $folder_location/data
bampaths_filename="$(basename -- $sample_bampaths)"

# Some of the scripts will expect tab-delimited metadata and species calls
sed -e "s/,/\t/g" $folder_location/data/$metadata_filename > $folder_location/data/sample_metadata_tabs.tsv

# We want to pull out the paths to the bamfile for each sample
bamfiles (){ 
	oldIFS=$IFS
	IFS=","
	read names
	while read $names; do
		echo -e ${sample_id}\\t${alignments_bam}
	done
	IFS=$oldIFS
}

cat $folder_location/data/$bampaths_filename | bamfiles > $folder_location/data/bampaths.csv


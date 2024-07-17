# This script prepares a folder in which a given sample set will be run
folder_location=$1
gsutil_path=$2
study_id=$3
sampleset=$4

# turn on extended globbing
shopt -s extglob

# Create the folder
mkdir -p $folder_location/data
# Download the metadata to the folder
for sset in ${sampleset[@]}; do
        sample_metadata=$gsutil_path/metadata/general/$sset/"samples.meta.csv
        gsutil cp $sample_metadata $folder_location/data/"$sset_samples.meta.csv"

        metadata_filename="$(basename -- $sample_metadata)"
        # Create a sample manifest from the metadata
        sed -e '1d' -e '2,$s/,.*//' $folder_location/data/$metadata_filename > > $folder_location/data/sample_manifest.txt
done

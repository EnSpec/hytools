#!/bin/bash

imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
hytools_dir=$(dirname ${imgspec_dir})

input="input"
mkdir -p output

# Get input paths for image correct
echo "Looking for input granule gzips and extracting if necessary..."
rfl_files=$(python ${imgspec_dir}/get_paths_from_granules.py -p rfl)
echo "Found input reflectance file(s): $rfl_files"
obs_ort_files=$(python ${imgspec_dir}/get_paths_from_granules.py -p obs_ort)
echo "Found input observation file(s): $obs_ort_files"

# Set ulimit according to ray recommendation
ulimit -n 8192

# Create image_correct_config
python ${imgspec_dir}/get_from_context.py image_correct_config > image_correct_config.json
echo "Created image_correct_config.json file from \"image_correct_config\" parameter"

if grep -q "file_type" image_correct_config.json; then
    python $imgspec_dir/update_config.py image_correct_config.json image_correct $rfl_files $obs_ort_files
    echo "Updated image_correct_config.json with input paths and num_cpus based on number of images"

    # Execute hytools image correction
    image_correct_cmd="python $hytools_dir/scripts/image_correct.py image_correct_config.json"
    echo "Executing cmd: $image_correct_cmd"
    ${image_correct_cmd}

    # Update rfl_files to refer to corrected output paths
    # Assume corrected output file names end with either "_topo" or "_brdf"
    rfl_files_arr=()
    for file in output/*{topo,brdf}; do
        if [[ $file != *\** ]]; then
            rfl_files_arr+=("$file")
        fi
    done
    rfl_files=$(printf "%s," "${rfl_files_arr[@]}" | cut -d "," -f 1-${#rfl_files_arr[@]})
    echo "Updated rfl_files based on output of image_correct step.  Found files:"
    echo $rfl_files

else
    echo "### image_correct_config.json is empty. Not running image correction step"
fi

# Create trait_estimate_config
python ${imgspec_dir}/get_from_context.py trait_estimate_config > trait_estimate_config.json
echo "Created trait_estimate_config.json file from \"trait_estimate_config\" parameter"

if grep -q "file_type" trait_estimate_config.json; then
    # Clone trait model repository
    trait_model_dir="trait_models"
    git clone -b $2 $1 $trait_model_dir

    python $imgspec_dir/update_config.py trait_estimate_config.json trait_estimate $rfl_files $obs_ort_files \
    $trait_model_dir
    echo "Updated trait_estimate_config.json with input paths, trait model paths, and num_cpus"

    # Execute hytools image correction
    trait_estimate_cmd="python $hytools_dir/scripts/trait_estimate.py trait_estimate_config.json"
    echo "Executing cmd: $trait_estimate_cmd"
    ${trait_estimate_cmd}
else
    echo "### trait_estimate_config.json is empty. Not running trait estimate step"
fi
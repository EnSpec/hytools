#!/bin/bash

imgspec_dir=$( cd "$(dirname "$0")" ; pwd -P )
hytools_dir=$(dirname $imgspec_dir)

input="input"
mkdir -p output

# Get input paths
rfl_files=$(python ${imgspec_dir}/get_paths_from_granules.py -p rfl)
echo "Found input reflectance file(s): $rfl_files"
obs_ort_files=$(python ${imgspec_dir}/get_paths_from_granules.py -p obs_ort)
echo "Found input observation file(s): $obs_ort_files"

# Create config.json
python get_from_context.py config_file > config.json
echo "Created config.json file from \"config_file\" parameter"
python $imgspec_dir/update_config.py config.json image_correct $rfl_files $obs_ort_files
echo "Updated config.json with input paths"

# Execute hytools image correction
image_correct_cmd="python $hytools_dir/scripts/image_correct.py config.json"
echo "Executing cmd: $image_correct_cmd"
${image_correct_cmd}
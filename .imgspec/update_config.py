#! /usr/bin/python

import glob
import json
import os
import sys


def main():

    # Read in config and input paths (assuming "envi" inputs for now)
    config_path = sys.argv[1]
    config_type = sys.argv[2]
    rfl_paths = sys.argv[3]
    obs_paths = sys.argv[4]

    with open(config_path, "r") as f:
        config_dict = json.loads(f.read())

    # Update input_files and anc_files fields
    config_dict['file_type'] = 'envi'
    aviris_anc_names = ['path_length','sensor_az','sensor_zn',
                        'solar_az', 'solar_zn','phase','slope',
                        'aspect', 'cosine_i','utc_time']
    images= [rfl for rfl in rfl_paths.split(",")]
    images.sort()
    config_dict["input_files"] = images

    config_dict["anc_files"] = {}
    anc_files = [obs for obs in obs_paths.split(",")]
    anc_files.sort()
    for i,image in enumerate(images):
        config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                    [[anc_files[i],a] for a in range(len(aviris_anc_names))]))

    # Update output_dir to local output directory
    output_dir = "output/"
    if config_type == "trait_estimate":
        config_dict["output_dir"] = output_dir
    else:
        config_dict["export"]["output_dir"] = output_dir

    # Update num_cpus based on number of images
    config_dict["num_cpus"] = len(images)

    if config_type == "trait_estimate":
        topo_coeffs_paths = sys.argv[5]
        brdf_coeffs_paths = sys.argv[6]

        topo_coeffs = [topo for topo in topo_coeffs_paths.split(",")]
        topo_coeffs.sort()

        brdf_coeffs = [brdf for brdf in brdf_coeffs_paths.split(",")]
        brdf_coeffs.sort()

        for corr in config_dict["corrections"]:
            config_dict[corr] = {}
            if corr == "topo":
                for i in range(len(topo_coeffs)):
                    config_dict[corr][images[i]] = topo_coeffs[i]
            if corr == "brdf":
                for i in range(len(brdf_coeffs)):
                    config_dict[corr][images[i]] = brdf_coeffs[i]

        imgspec_dir = os.path.dirname(os.path.abspath(__file__))
        models = glob.glob(os.path.join(imgspec_dir, "trait_models", "*json"))
        config_dict["trait_models"] = models

    with open(config_path, "w") as f:
        f.write(json.dumps(config_dict))


if __name__ == "__main__":
    main()
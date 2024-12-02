## FlexBRDF without Ray

In some distributed system like [High Throughput Computing system](https://chtc.cs.wisc.edu/), the workers are not sharing memory, or computation or storage resource are limited. In such situations, not enough computer cores can be assigned to the FlexBRDF process that matches the total number of flightlines in the same group. Ray cannot be run easily with limited umber of CPU cores.

To solve that in the FlexBRDF context, samples that needed for BRDF modeling from each flightline can be extracted and saved individually, and then transferred to the same storage location later for the final step of BRDF coefficient estimation. The extraction step can be either executed sequentially or parallelly, depending on the system.

Extracted samples are stored in HDF files in the current solution, which is one h5 file for each flightline in the group. Once all extractions are done, all h5 files are gathered and used for the final step.

Take HTC as an example, DAG management can be used to monitor and control the work flow.

Two scripts will be used consecutively for this purpose
([image_correct_get_sample_chtc.py](../scripts/no_ray/image_correct_get_sample_chtc.py) and [image_correct_combine_sample_chtc.py](../scripts/no_ray/image_correct_combine_sample_chtc.py)).

### Extraction

Target flightline can be set with the last commandline parameter. The order is determined by the order of the ```input_files``` in the configuration json file. In this step, each single run can be independent. The output h5 file will be stored in the output path defined in the configuration file, which may or may not necessarily be the same for each run. If ```topo``` is enabled, the resultant coefficient json file will also be stored there, and the extracted reflectance samples in h5 will be topographically corrected. 

Please make sure the bands for computing NDVI are in the good band range of the configuration file.

```bash
# first line
python ./scripts/no_ray/image_correct_get_sample_chtc.py path/to/the/configuration/json/file  0 
# second line
python ./scripts/no_ray/image_correct_get_sample_chtc.py path/to/the/configuration/json/file  1 
# third line
python ./scripts/no_ray/image_correct_get_sample_chtc.py path/to/the/configuration/json/file  2 
# so forth
... ...
```

### Combination

After all the independent jobs are finished, the final step is to combine all samples and continue the BRDF modeling. 

```bash
python ./scripts/no_ray/image_correct_combine_sample_chtc.py path/to/the/configuration/json/file folder/of/the/h5files/in/the/same/group
```

All brdf coefficient json files will be stored in the same path specified in the configuration json file.

### Script to run the two steps together

If scripts are not run in a HTC-like system, they can still be run in a single machine with multiple CPU cores. This script is a simplified version of the workflow for combining two steps of FlexBRDF ([run_single_process_merge.py](../scripts/no_ray/run_single_process_merge.py)). It requires "*path/to/the/configuration/json/file*" and the total number of flightlines in the group as inputs. It initiates multiple instances of the script [image_correct_get_sample_chtc.py](../scripts/no_ray/image_correct_get_sample_chtc.py)

```python
import sys, os
import multiprocessing
import subprocess, json

exec_str="python ./script/no_ray/image_correct_get_sample_chtc.py "

merge_str="python ./script/no_ray/image_correct_combine_sample_chtc.py {} {}"

def run_command(command):
    print(command)
    subprocess.run(command,shell=True) 


def main():

    config_file = sys.argv[1]
    total_count = int(sys.argv[2])
    worker_count = min(os.cpu_count()-1,total_count)
    
    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)
        h5_folder=config_dict["export"]["output_dir"]

        pool = multiprocessing.Pool(processes=worker_count)

        commands =  [f"{exec_str} {config_file} {order}" for order in range(total_count)]
        pool.map(run_command, commands)
        pool.close()
        pool.join()  # Wait for all subprocesses to finish

        print('All extractions are done.')
        
        # Final step to use all pixels in the group to estimate BRDF coefficients
        subprocess.run(merge_str.format(config_file,h5_folder),shell=True) 

if __name__== "__main__":
    main()
```

### Script to run the three steps together for Subgrouped-TOPO correction 
There is also a script for combining three steps of FlexBRDF with subgroup for TOPO ([run_single_process_topogroup_merge_v3.py](../scripts/no_ray/run_single_process_topogroup_merge_v3.py)) It controls and monitors the workflow to run three related scripts ([image_correct_get_raw_sample_chtc.py](../scripts/no_ray/image_correct_get_raw_sample_chtc.py) for extracting raw samples, [image_correct_combine_topo_sample_chtc.py](../scripts/no_ray/image_correct_combine_topo_sample_chtc.py) for gathering samples in the same TOPO subgroup, and the final script for merging them for BRDF [image_correct_combine_sample_chtc.py](../scripts/no_ray/image_correct_combine_sample_chtc.py)). 


### Export corrected images

After the previous two steps for TOPO/BRDF model estimation, user can export corrected images with those precomputed coefficients. 

This script is a simplified version of the workflow for exporting images and masks ([run_single_process_export.py](../scripts/no_ray/run_single_process_export.py)). It requires "*path/to/the/configuration/json/file*" and the total number of flightlines in the group as inputs. It initiates multiple instances of the script [image_correct_export_image.py](../scripts/no_ray/image_correct_export_image.py)


### Export trait estimation with corrected image

The correction coefficients can also be used for trait mapping. If the any of the three corrections is enabled, trait prediction will be applied to corrected pixels.

This script is a simplified version of the workflow for exporting images and masks ([run_single_process_trait.py](../scripts/no_ray/run_single_process_trait.py)). It requires "*path/to/the/trait/export/configuration/json/file*", the total number of flightlines in the group, and the total number of traits as inputs. It initiates multiple instances of the script [trait_estimate_inde.py](../scripts/no_ray/trait_estimate_inde.py)


import sys, os, time
import multiprocessing
import subprocess, json

exec_str="python ../noray/image_correct_get_raw_sample_chtc.py "

topo_group_str = "python ../noray/image_correct_combine_topo_sample_chtc.py {} {} {}"

merge_str="python ../noray/image_correct_combine_sample_chtc.py {} {}"

def parse_group_info(group_dict,full_img_list):
    out_group_dict={}
    out_group_dict_short={}
    for img_name in full_img_list:
        subgroup_name = group_dict[img_name]
        if subgroup_name in out_group_dict:
            out_group_dict[subgroup_name]+=[img_name]
            out_group_dict_short[subgroup_name]+=[full_img_list.index(img_name)] 
        else:
            out_group_dict[subgroup_name]=[img_name]
            out_group_dict_short[subgroup_name]=[full_img_list.index(img_name)]

    return out_group_dict, out_group_dict_short

def run_command(command):
    subprocess.run(command,shell=True)


def run_step2(params):
    config_file,h5_folder,current_group_name=params
    #print(topo_group_str.format(config_file,h5_folder,current_group_name))
    subprocess.run(topo_group_str.format(config_file,h5_folder,current_group_name),shell=True)


def main():

    config_file = sys.argv[1]
    h5_folder =  sys.argv[2]

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

        image_list = config_dict["input_files"]
        total_count = int(len(image_list))

        subgroup=config_dict['topo']['subgroup']
        group_meta_dict, worker_unfinished=parse_group_info(subgroup,image_list)

        level1_worker_count = min(os.cpu_count()-1,total_count)

        print(level1_worker_count)

        with multiprocessing.Pool(processes=level1_worker_count) as  pool:

            subgroup_status_dict={}
            for sub_group_name in group_meta_dict:
                step1_status=[]
                for order_in_full_list in worker_unfinished[sub_group_name]:
                    command = f"{exec_str} {config_file} {order_in_full_list} 1"
                    step1_status += [pool.apply_async(run_command, args=(command,))]
                subgroup_status_dict[sub_group_name] = step1_status

            level1_total_unfinished = total_count
            while level1_total_unfinished:
                # check for finished subgroup every 10 sec. If the raw pixel extraction is done for a subgroup, start step2-group TOPO. 
                temp_unfinished_count=0
                subgroup_list = list(subgroup_status_dict.keys())
                for gp_name in subgroup_list:
                    sub_count=sum([not r.ready() for r in subgroup_status_dict[gp_name]])

                    if sub_count==0:
                        pool.apply_async(run_step2, args=((config_file,h5_folder,gp_name),))
                        del subgroup_status_dict[gp_name]
                    temp_unfinished_count += sub_count
                level1_total_unfinished=temp_unfinished_count

                #print(f"=={level1_total_unfinished}/{total_count} tasks Level 1 pretopo")
                time.sleep(10.0)

        print('All extraction is done.')
        subprocess.run(merge_str.format(config_file,h5_folder),shell=True) 

if __name__== "__main__":
    main()

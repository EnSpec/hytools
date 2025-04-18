
import sys, os
import multiprocessing
import subprocess, json

exec_str="python ../no_ray/image_correct_get_sample_chtc.py "

merge_str="python ../no_ray/image_correct_combine_sample_chtc.py {} {}"

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

        if total_count > len(config_dict["input_files"]):
            print("Out of upper bound")
            return

        pool = multiprocessing.Pool(processes=worker_count)

        commands =  [f"{exec_str} {config_file} {order}" for order in range(total_count)]
        pool.map(run_command, commands)
        pool.close()
        pool.join()  # Wait for all subprocesses to finish

        print('All extractions are done.')

        subprocess.run(merge_str.format(config_file,h5_folder),shell=True) 

if __name__== "__main__":
    main()
    
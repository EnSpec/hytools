
import sys, os 
import multiprocessing
import subprocess, json

exec_str="python /home/ye6/hys_test/gdrive_test/trait_estimate_inde.py "   

def run_command(command):
    print(command)
    subprocess.run(command,shell=True) 

def main():

    config_file = sys.argv[1]

    total_img_count = int(sys.argv[2])
    total_trait_count = int(sys.argv[3])

    total_count = total_img_count*total_trait_count
    worker_count = min(os.cpu_count()-1,total_count)

    with open(config_file, 'r') as outfile:
        config_dict = json.load(outfile)

        pool = multiprocessing.Pool(processes=worker_count)

        if (total_img_count > len(config_dict["input_files"])) or (total_trait_count> len(config_dict["input_files"])):
            print("Out of upper bound")
            return

        param_list = []
        for img_i in range(total_img_count):
            for trait_j in range(total_trait_count):
                param_list+=[[img_i,trait_j]]

        commands =  [f"{exec_str} {config_file} {param_order[0]} {param_order[1]}" for param_order in param_list]  #[f"{exec_str} 0", f"{exec_str} 1"]
        pool.map(run_command, commands)
        pool.close()
        pool.join()  # Wait for all subprocesses to finish

        print('All traits are done.')


if __name__== "__main__":
    main()
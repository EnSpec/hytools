
import sys, os 
import multiprocessing
import subprocess, json

exec_str="python ~/hys_test/gdrive_test/image_correct_export_image.py "   

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

        commands =  [f"{exec_str} {config_file} {order}" for order in range(total_count)]  #[f"{exec_str} 0", f"{exec_str} 1"]
        pool.map(run_command, commands)
        pool.close()
        pool.join()  # Wait for all subprocesses to finish

        print('All image export is done.')

if __name__== "__main__":
    main()
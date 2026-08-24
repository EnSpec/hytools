
import os
import sys
import tarfile

intar = sys.argv[1]

indir= sys.argv[2]

out_dir = sys.argv[3]

AVNG_ftp_path = "https://avng.jpl.nasa.gov/avng/y{}/"
AVC_ftp_path = "https://popo.jpl.nasa.gov/avcl/y{}_data/"



def generate_download_tag(intar):

    year_tag, ftp_path, file_tag_rfl, file_tag_obs = None, None, None, None

    if intar.startswith('ang'):
        # AVNG
        year_tag = intar[5:7]
        #print(year_tag)

        if year_tag>='15':
            ftp_path = AVNG_ftp_path.format(year_tag)
        else:
            print("Not the right year.")
            sys.exit(1)

        file_tag_rfl = f"{intar}rfl.tar.gz"
        file_tag_obs =f"{intar}.tar.gz"

    elif intar.startswith('f'):
        # AVC
        year_tag = intar[1:3]
        #print(year_tag)
        if year_tag>'11':
            ftp_path = AVC_ftp_path.format(year_tag)

        else:
            print("Not the right year. Check local...")
            #sys.exit(1)  

        #file_tag_rfl = "{}rfl.tar.gz".format(intar)
        file_tag_rfl = f"{intar}_refl.tar.gz"
        file_tag_obs = f"{intar}.tar.gz"

    else:
        sys.exit(1)

    return year_tag, ftp_path, file_tag_rfl, file_tag_obs


def extract_refl(intar, indir, out_dir, file_tag_rfl):
    # L2 reflectance
    #tar = tarfile.open(indir+'/'+intar+'rfl.tar.gz', "r:gz")
    #tar = tarfile.open(indir+'/'+file_tag_rfl, "r:gz")
    with tarfile.open(indir+'/'+file_tag_rfl, "r:gz") as tar:
        tar_folder_name=None
        version_tag=''
        for tarinfo in tar:
            print(tarinfo.name,tarinfo.isdir(),tarinfo.isreg())
            #print(tarinfo.name, "is", tarinfo.size, "bytes in size and is", end="")
            if intar.startswith('ang'):
                if tarinfo.isreg():
                    dir_name_list = tarinfo.name.split('_')
                    if dir_name_list[-1].startswith('v'):
                        version_tag = dir_name_list[-1]
                    elif dir_name_list[-2].startswith('v'):
                        version_tag = dir_name_list[-2]

                    #print(" a regular file.")
                    if tarinfo.name.endswith('corr_'+version_tag+'_img') or tarinfo.name.endswith('rfl_'+version_tag+'_img'):
                        print(tarinfo.name)

                        #print('_'.join(tarinfo.name.split('_')[-3:]))
                        base_name = intar+'_'+'_'.join(tarinfo.name.split('_')[-3:])
                        print(base_name)

                        tar.extract(tarinfo,path=out_dir)
                        tar.extract(tarinfo.name+'.hdr',path=out_dir)

                        os.rename(os.path.join(out_dir,tarinfo.name+'.hdr'),os.path.join(out_dir,base_name+'.hdr'))
                        os.rename(os.path.join(out_dir,tarinfo.name),os.path.join(out_dir,base_name))

                        if tar_folder_name is not None:
                            os.rmdir(os.path.join(out_dir, tar_folder_name))

                        break

                elif tarinfo.isdir():
                    print(" a directory.")
                    #print(tarinfo.name.split('_'))
                    tar_folder_name = tarinfo.name
                    dir_name_list = tarinfo.name.split('_')
                    if dir_name_list[-1].startswith('v'):
                        version_tag = dir_name_list[-1]
                    elif dir_name_list[-2].startswith('v'):
                        version_tag = dir_name_list[-2]
                    #break
                else:
                    print(" something else.")
                #print(version_tag)
              
            elif intar.startswith('f'):
                if tarinfo.isreg():
                    #print(" a regular file.")
                    #(tarinfo.name.split(intar)[-1])
                    refl_tag = tarinfo.name.split(intar)[-1]
                    #print(refl_tag)
                    #continue
                    if refl_tag.startswith('_corr_') or refl_tag.startswith('_rfl_') or refl_tag.startswith('rdn_refl_') or refl_tag.startswith('rdn_rfl_'):
                        #print(tarinfo.name)
                        #print('_'.join(tarinfo.name.split('_')[-3:]))
                        base_name = (intar+refl_tag).split('.hdr')[0]
                        print(base_name)

                        tar.extract(tarinfo.name.split('.hdr')[0],path=out_dir)
                        tar.extract(tarinfo.name.split('.hdr')[0]+'.hdr',path=out_dir)
                        
                        os.rename(os.path.join(out_dir,tarinfo.name.split('.hdr')[0]+'.hdr'),os.path.join(out_dir,base_name+'.hdr'))
                        os.rename(os.path.join(out_dir,tarinfo.name.split('.hdr')[0]),os.path.join(out_dir,base_name))

                        os.rmdir(os.path.join(out_dir, tar_folder_name))

                        break

                elif tarinfo.isdir():
                    #print(" a directory.")
                    #print(tarinfo.name.split('_'))
                    tar_folder_name = tarinfo.name
                    version_tag = None
                else:
                    print(" something else.")
                    version_tag = None
        #tar.close()
        #return base_name

def extract_obs_ort(intar, indir, out_dir, file_tag_obs):
# L1 obs_ort
#tar = tarfile.open(indir+'/'+intar+'.tar.gz', "r:gz")
    with tarfile.open(indir+'/'+file_tag_obs, "r:gz") as tar:
    #tar = tarfile.open(indir+'/'+file_tag_obs, "r:gz")
        tar_folder_name = None
        version_tag = ''
        for tarinfo in tar:
            #print(tarinfo.name)
            #print(tarinfo.name, "is", tarinfo.size, "bytes in size and is", end="")
            if tarinfo.isreg():
                #print(" a regular file.")
                #if tarinfo.name.endswith(version_tag+'_obs_ort') :
                if tarinfo.name.endswith('obs_ort') or tarinfo.name.endswith('OBS_ORT'):
                    print(tarinfo.name)
                    #print('_'.join(tarinfo.name.split('_')[-4:]))
                    #base_name = intar+'_'+'_'.join(tarinfo.name.split('_')[-4:])
                    base_name = intar+version_tag+'_obs_ort' #tarinfo.name.split(os.path.sep)[-1]

                    tar.extract(tarinfo,path=out_dir)
                    tar.extract(tarinfo.name+'.hdr',path=out_dir)

                    os.rename(os.path.join(out_dir,tarinfo.name+'.hdr'),os.path.join(out_dir,base_name+'.hdr'))
                    os.rename(os.path.join(out_dir,tarinfo.name),os.path.join(out_dir,base_name))

                    if tar_folder_name is not None:
                        os.rmdir(os.path.join(out_dir, tar_folder_name))

                    break
            elif tarinfo.isdir():
                #print(" a directory.")
                print(tarinfo.name.split('_'))
                print(tarinfo.name.split(intar)[-1])
                version_tag=tarinfo.name.split(intar)[-1]

                tar_folder_name = tarinfo.name

                #break
            else:
                print(" something else.")
                #version_tag = ''
                #tar_folder_name = ''
        #print(version_tag)
        print(base_name)
      #tar.close()
    return base_name


year_tag, ftp_path, file_tag_rfl, file_tag_obs = generate_download_tag(intar)
print(year_tag, ftp_path, file_tag_rfl, file_tag_obs)

extract_refl(intar, indir, out_dir, file_tag_rfl)
extract_obs_ort(intar, indir, out_dir, file_tag_obs)

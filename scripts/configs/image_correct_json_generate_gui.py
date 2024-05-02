
# python c:/mydata/image_correct_json_generate_gui.py

import os
import json
import glob
import numpy as np

import tkinter as tk
from tkinter.filedialog import asksaveasfilename, askdirectory


def fill_config(images, anc_files,out_coef_dir,img_file_type,corr_list, flag_pre_compute=False, topo_coeff = [], brdf_coeff=[]):

  config_dict = {}

  #Only coefficients for good bands will be calculated
  config_dict['bad_bands'] =[[300,400],[1337,1430],[1800,1960],[2450,2600]]
  
  config_dict['file_type'] = img_file_type #'envi'
  #aviris_anc_names = ['path_length','sensor_az','sensor_zn',
  #                    'solar_az', 'solar_zn','phase','slope',
  #                    'aspect', 'cosine_i','utc_time']
  aviris_anc_names = ['sensor_az','sensor_zn',
                      'solar_az', 'solar_zn']                    
  images.sort()
  config_dict["input_files"] = images

  if img_file_type=='envi':
    config_dict["anc_files"] = {}
    anc_files.sort()
    for i,image in enumerate(images):
        config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                    [[anc_files[i],a] for a in range(len(aviris_anc_names))]))  
  
  config_dict['export'] = {}
  config_dict["topo"] =  {}
  config_dict["brdf"]  = {}
  
  if flag_pre_compute:
    config_dict['export']['coeffs']  = False
    config_dict['export']['image']  = True    
  else:  
    config_dict['export']['coeffs']  = True
    config_dict['export']['image']  = False
  
  config_dict['export']['masks']  = False
  config_dict['export']['subset_waves']  =  [660,550,440] #[440,550,660,850] #
  config_dict['export']['output_dir'] = out_coef_dir
  print('_'.join(corr_list))
  config_dict['export']["suffix"] = '_'.join(corr_list)  # 'brdf'  
  
  config_dict["corrections"]  =  corr_list # ['brdf'] 
  
  
  if flag_pre_compute and len(topo_coeff)==len(images):
    
    config_dict["topo"]['type'] =  'precomputed'
    topo_files = sorted(topo_coeff)
    #print(dict(zip(images, topo_files)))
    config_dict["topo"]['coeff_files'] = dict(zip(images, topo_files))
  else:  
    config_dict["topo"]['type'] =  'scs+c'
    
    config_dict["topo"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                               'min': 0.05,'max': 1.0}]
                                     ]
    config_dict["topo"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                               'min': 0.05,'max': 1.0}]]
    config_dict["topo"]['c_fit_type'] =   'nnls' #'ols' #'nnls' #
  
  if flag_pre_compute and len(brdf_coeff)==len(images):  
    
    config_dict["brdf"]['type'] =  'precomputed'
    brdf_files = sorted(brdf_coeff)
    config_dict["brdf"]['coeff_files'] = dict(zip(images, brdf_files))
  else:    
    # Options are 'line','scene', or a float for a custom solar zn
    # Custom solar zenith angle should be in radians
    config_dict["brdf"]['solar_zn_type'] ='scene'

    #----------------------
    # ## Flex BRDF configs
    # ##------------------
    config_dict["brdf"]['type'] =  'flex'
    config_dict["brdf"]['grouped'] =  True
    config_dict["brdf"]['geometric'] = 'li_sparse_r'
    config_dict["brdf"]['volume'] = 'ross_thick'
    config_dict["brdf"]["b/r"] = 2.5
    config_dict["brdf"]["h/b"] = 2
    config_dict["brdf"]['sample_perc'] = 0.1
    config_dict["brdf"]['interp_kind'] = 'linear'
    config_dict["brdf"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                 'min': 0.05,'max': 1.0}],
                                        ['kernel_finite',{}],
                                        ['ancillary',{'name':'sensor_zn',
                                                      'min':np.radians(2),'max':'inf' }]
                                       ]


    config_dict["brdf"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                  'min': 0.05,'max': 1.0}]]

    # ## Flex dynamic NDVI params
    config_dict["brdf"]['bin_type'] = 'dynamic'
    config_dict["brdf"]['num_bins'] = 18
    config_dict["brdf"]['ndvi_bin_min'] = 0.05
    config_dict["brdf"]['ndvi_bin_max'] = 1.0
    config_dict["brdf"]['ndvi_perc_min'] = 10
    config_dict["brdf"]['ndvi_perc_max'] = 95
  

  config_dict["resample"]  = False
  config_dict['num_cpus'] = len(images)  

  return config_dict

'''  
def update_corr_list(corr_list_):
   #print(chk_topo.get(),chk_brdf.get())
   corr_list_ = ['topo']*chk_topo.get()+['brdf']*chk_brdf.get()
   #print(corr_list)
'''  
  
def gen_config(entry_outdir, entry_outjson,img_list_out, obs_list_out, radio_f_type, corr_list,chk_precompute, topo_list_out,brdf_list_out):  

  outdir_name = entry_outdir.get()
  out_json = entry_outjson.get()
  images = img_list_out['text'].split('\n')
  anc_files = obs_list_out['text'].split('\n')
  img_file_type = str(radio_f_type.get()).lower()
  corr_list_str = ['topo']*(corr_list[0].get()) + ['brdf']*(corr_list[1].get())
  flag_pre_compute = bool(chk_precompute.get())
  #print(flag_pre_compute,chk_precompute.get())
  if flag_pre_compute:
    topo_json_list = topo_list_out['text'].split('\n')
    brdf_json_list = brdf_list_out['text'].split('\n')

    return_json_dict = fill_config(images, anc_files,outdir_name,img_file_type,corr_list_str,flag_pre_compute=flag_pre_compute, topo_coeff = topo_json_list, brdf_coeff=brdf_json_list)
  else:
  
    return_json_dict = fill_config(images, anc_files,outdir_name,img_file_type,corr_list_str)
  
  with open(out_json, 'w') as outfile:
    json.dump(return_json_dict,outfile,indent=3)
    window.title(f"File saved- {out_json}")
  
def save_file(txt_out_json):
    """Save the current file as a new file."""
    filepath = asksaveasfilename(
        defaultextension=".json",
        filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")],
    )

    if not filepath:
        return

    txt_out_json.delete(0,tk.END)
    txt_out_json.insert(0,filepath)

    window.title(f"Export Configuration JSON  - {filepath}")

def open_folder(out_component):
    """Open a file for editing."""
    in_img_dir = askdirectory()
    if not in_img_dir:
        return
     
    #print(in_img_dir)     
    out_component.delete(0, tk.END)    
    out_component.insert(0,in_img_dir)    

    
    window.title(f"Folder- {in_img_dir}")    
    
def open_file(out_component, list_component, radio_f_type=None, pattern=None):
    """Open a file for editing."""
    in_img_dir = askdirectory()
    if not in_img_dir:
        return
        
    #out_component.delete("1.0", tk.END)    
    #out_component.insert("1.0",in_img_dir)    
    out_component["text"] = in_img_dir
    if pattern is None:
      img_file_type = str(radio_f_type.get()).lower()
      file_ext_list = {'envi':'*img','neon':'*.h5'}
      pattern = file_ext_list[img_file_type]

      
    #if pattern is None:
    #  return  in_img_dir   
        
    file_list = sorted(glob.glob(in_img_dir+'/'+pattern))

    if len(file_list)==0:
      list_component['text']= 'No files selected'
      return 
        
    list_component['text']= '\n'.join([ os.path.normpath(x) for x in file_list])
    #print(list_component['text'])
    
    window.title(f"Folder- {in_img_dir}")

window = tk.Tk()
window.title("Setup image correction configuration file")

#window.rowconfigure(0, minsize=400, weight=1)
#window.columnconfigure(1, minsize=600, weight=1)


#txt_edit = tk.Text(window)

frm_img_buttons = tk.Frame(window, relief=tk.RAISED, bd=2)
frm_img_buttons.columnconfigure(1, minsize=300, weight=1)
frm_img_buttons.rowconfigure(1, minsize=50, weight=1)

frm_obs_buttons = tk.Frame(window, relief=tk.RAISED, bd=2)
frm_obs_buttons.columnconfigure(1, minsize=300, weight=1)
frm_obs_buttons.rowconfigure(1, minsize=50, weight=1)

frm_out_buttons = tk.Frame(window, relief=tk.RAISED, bd=2)
frm_out_buttons.columnconfigure(1, weight=1) #, minsize=300

frm_out_json_buttons = tk.Frame(window, relief=tk.RAISED, bd=2)
frm_out_json_buttons.columnconfigure(1,  weight=1) #minsize=300,

frm_file_type = tk.Frame(window, relief=tk.RAISED, bd=2)
frm_corr_type = tk.Frame(window, relief=tk.RAISED, bd=2)

frm_precomp = tk.Frame(window, relief=tk.RAISED, bd=2,highlightbackground="grey",highlightthickness=5,padx=5, pady=5)
frm_precomp.columnconfigure(1,weight=1) 
frm_precomp.rowconfigure(1, minsize=50, weight=1) 
frm_pre_topo_buttons = tk.Frame(frm_precomp, relief=tk.RAISED, bd=2)
frm_pre_topo_buttons.columnconfigure(1, minsize=280, weight=1)
#frm_pre_topo_buttons.rowconfigure(1, minsize=50, weight=1)

frm_pre_brdf_buttons = tk.Frame(frm_precomp, relief=tk.RAISED, bd=2)
frm_pre_brdf_buttons.columnconfigure(1, minsize=280, weight=1)
#frm_pre_brdf_buttons.rowconfigure(1, minsize=50, weight=1)

frm_final_gen = tk.Frame(window, relief=tk.RAISED, bd=2,highlightbackground="grey",highlightthickness=5)
frm_final_gen.rowconfigure(0, minsize=50, weight=0) 
frm_final_gen.columnconfigure(0,  minsize=100, weight=1) #
frm_final_gen.columnconfigure(1,  minsize=50, weight=0)
frm_final_gen.columnconfigure(2,  minsize=100, weight=1)

label_img_dir = tk.Label(frm_img_buttons,text="image", fg="white", bg="black")
label_obs_dir = tk.Label(frm_obs_buttons,text="image", fg="white", bg="black")
txt_outdir = tk.Entry(frm_out_buttons)

img_list_out = tk.Label(frm_img_buttons,bg="grey",anchor="w")
obs_list_out = tk.Label(frm_obs_buttons,bg="grey",anchor="w")


btn_open1 = tk.Button(frm_img_buttons, text="Image Folder...", command=lambda: open_file(label_img_dir,img_list_out,radio_f_type=f_type) )  #'*img'
btn_open2 = tk.Button(frm_obs_buttons, text="Obs_ort Folder...", command=lambda: open_file(label_obs_dir,obs_list_out,pattern='*obs_ort'))
btn_open3 = tk.Button(frm_out_buttons, text="Output coeff Folder...", command=lambda: open_folder(txt_outdir)) #, command=open_file(txt_edit, None)
#btn_open1 = tk.Button(frm_img_buttons, text="Image Folder...", command=lambda: open_file(label_img_dir,img_list_out,radio_f_type=f_type) )  #'*img'
#btn_open2 = tk.Button(frm_obs_buttons, text="Obs_ort Folder...", command=lambda: open_file(label_obs_dir,obs_list_out,pattern='*obs_ort'))

btn_open_topo = tk.Button(frm_pre_topo_buttons, text="TOPO json Folder...", command=lambda: open_file(label_pre_topo_dir,topo_list_out,pattern='*topo_coeffs*.json'))
btn_open_brdf = tk.Button(frm_pre_brdf_buttons, text="BRDF json Folder...", command=lambda: open_file(label_pre_brdf_dir,brdf_list_out,pattern='*brdf_coeffs*.json'))

txt_out_json = tk.Entry(frm_out_json_buttons)

btn_save = tk.Button(frm_out_json_buttons, text="Save As...", command=lambda: save_file(txt_out_json)) 

btn_open1.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
label_img_dir.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
img_list_out.grid(row=1, columnspan=2, sticky="nsew", padx=5, pady=5)

btn_open2.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
label_obs_dir.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
obs_list_out.grid(row=1, columnspan=2, sticky="nsew", padx=5, pady=5)

btn_open3.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
txt_outdir.grid(row=0, column=1, sticky="ew",  padx=5, pady=5)
btn_save.grid(row=0, column=0, sticky="ew", padx=5)
txt_out_json.grid(row=0, column=1, sticky="ew",  padx=5, pady=5)

f_type = tk.StringVar(value="envi")
ev_btn = tk.Radiobutton(frm_file_type, text='ENVI (*img)', variable=f_type, value='envi')
h5_btn = tk.Radiobutton(frm_file_type, text='NEON HDF5 (*.h5)', variable=f_type, value='neon')
ev_btn.grid(row=0, column=0, sticky="ew", padx=5)
h5_btn.grid(row=0, column=1, sticky="ew", padx=5)

chk_topo = tk.IntVar(value=0)
chk_brdf = tk.IntVar(value=1)

chk_topo_btn = tk.Checkbutton(frm_corr_type, text='TOPO', variable=chk_topo, onvalue=1, offvalue=0)  #, command=lambda: update_corr_list(corr_list))
chk_brdf_btn = tk.Checkbutton(frm_corr_type, text='BRDF', variable=chk_brdf, onvalue=1, offvalue=0) #, command=lambda: update_corr_list(corr_list))
chk_topo_btn.grid(row=0, column=0, sticky="ew", padx=5)
chk_brdf_btn.grid(row=0, column=1, sticky="ew", padx=5)
corr_list = [chk_topo,chk_brdf]

chk_precompute = tk.IntVar(value=0)
coeff_list= []
label_precomput = tk.Checkbutton(frm_precomp, text='Load Precomputed Coefficients', variable=chk_precompute, onvalue=1, offvalue=0,anchor="center")  #tk.Label(frm_precomp,text="Precomputed Coefficients")   
label_precomput.grid(row=0, columnspan=2, sticky="ew",  padx=2, pady=2)
btn_open_topo.grid(row=0, column=0, sticky="ew", padx=2, pady=2)
btn_open_brdf.grid(row=0, column=0, sticky="ew", padx=2, pady=2)
label_pre_topo_dir = tk.Label(frm_pre_topo_buttons,text="topo json", fg="white", bg="black")
label_pre_topo_dir.grid(row=0, column=1, sticky="ew",  padx=2, pady=2)
label_pre_brdf_dir = tk.Label(frm_pre_brdf_buttons,text="brdf json", fg="white", bg="black")
label_pre_brdf_dir.grid(row=0, column=1, sticky="ew",  padx=2, pady=2)
topo_list_out = tk.Label(frm_pre_topo_buttons,bg="grey",anchor="w")
topo_list_out.grid(row=1, columnspan=2, sticky="ew",  padx=2,pady=2)
brdf_list_out = tk.Label(frm_pre_brdf_buttons,bg="grey",anchor="w")
brdf_list_out.grid(row=1, columnspan=2, sticky="ew",  padx=2,pady=2)
frm_pre_topo_buttons.grid(row=1, column=0, sticky="ew")
frm_pre_brdf_buttons.grid(row=1, column=1, sticky="ew")


#img_list_out = tk.Label(frm_img_buttons,bg="grey",anchor="w")
#obs_list_out = tk.Label(frm_obs_buttons,bg="grey",anchor="w")



btn_gen = tk.Button(frm_final_gen, text="Generate", font=("Calibri",12,"bold"), command=lambda: gen_config(txt_outdir,txt_out_json,img_list_out, obs_list_out,f_type,corr_list,chk_precompute,topo_list_out, brdf_list_out))
#btn_gen.place(relx=.5, rely=.5,anchor= 'e')
btn_gen.grid(row=0,column=1,sticky="wens") #anchor='center', 
#btn_gen.place(relx=0.5, rely=0.95, anchor=tk.CENTER)

frm_img_buttons.grid(row=0, column=0, sticky="ew")
frm_obs_buttons.grid(row=0, column=1, sticky="ew")
frm_out_buttons.grid(row=1, columnspan=2, sticky="we")
frm_file_type.grid(row=2,column=0, sticky="ew")
frm_corr_type.grid(row=2,column=1, sticky="ew")
frm_precomp.grid(row=3,columnspan=2, sticky="nsew")
frm_out_json_buttons.grid(row=4, columnspan=2, sticky="ew")
frm_final_gen.grid(row=6, columnspan=2, sticky="nsew") #, sticky="nsew"

window.mainloop()
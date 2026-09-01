import tkinter as tk
from tkinter import ttk, filedialog, messagebox
#from tkinter import ttk, filedialog, messagebox, scrolledtext
import json
import os
import csv
from datetime import datetime, timezone
import numpy as np

def generate_config(setting_dict):


    feedback_message = []
    corr_list = setting_dict["processing_pipeline"]["corrections"]
    ascii_loader_bool = setting_dict["input_datasets"]["ascii_loader_enabled"]

    config_dict = {}
    config_dict["job_metadata"] = setting_dict["job_metadata"]
    config_dict['bad_bands'] = setting_dict["processing_pipeline"]["bad_bands"]
    config_dict["file_type"] = setting_dict["processing_pipeline"]["data_type"]

    if ascii_loader_bool: # do not sort, use ascii order
        images = setting_dict["input_datasets"]["refl_images"]
    else:
        images = sorted(setting_dict["input_datasets"]["refl_images"])

    config_dict["input_files"] = images
    config_dict['num_cpus'] = len(images)

    if len(images)<1:
        feedback_message+=['// No reflectance images is selected.']

    aviris_anc_names = ['path_length','sensor_az','sensor_zn',
                        'solar_az', 'solar_zn','phase','slope',
                        'aspect', 'cosine_i','utc_time']

    if config_dict["file_type"] in ["envi","emit"]:

        config_dict["anc_files"] = {}

        if ascii_loader_bool: # do not sort, use ascii order
            anc_files = setting_dict["input_datasets"]["anc_files"]
        else:    
            anc_files = sorted(setting_dict["input_datasets"]["anc_files"])

        if 'topo' in corr_list or 'brdf' in corr_list:
            if not len(images)==len(anc_files):
                feedback_message+=['// Reflectance images and ancillary images do not match in number.']
            elif len(anc_files)==0:
                feedback_message+=['// Ancillary images required for TOPO or BRDF correction, but none is provided.']
            else:
                for i,image in enumerate(images):
                    config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                                [[anc_files[i],a] for a in range(len(aviris_anc_names))]))

    elif config_dict["file_type"] =="neon":
        pass
    elif config_dict["file_type"] in ["ncav"]:

        ncav_anc_names = ['path_length','to_sensor_azimuth','to_sensor_zenith',
                        'to_sun_azimuth', 'to_sun_zenith','solar_phase','slope',
                        'aspect', 'cosine_i','utc_time']

        config_dict["anc_files"] = {}

        if ascii_loader_bool: # do not sort, use ascii order
            anc_files = setting_dict["input_datasets"]["anc_files"]
        else:
            anc_files = sorted(setting_dict["input_datasets"]["anc_files"])

        if 'topo' in corr_list or 'brdf' in corr_list:
            if not len(images)==len(anc_files):
                feedback_message+=['// Reflectance images and ancillary images do not match in number.']
            elif len(anc_files)==0:
                feedback_message+=['// Ancillary images required for TOPO or BRDF correction, but none is provided.']
            else:
                for i,image in enumerate(images):
                    config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                                [[anc_files[i],ncav_anc_names[a]] for a in range(len(aviris_anc_names))]))

    elif config_dict["file_type"] in ["tanager"]:

        external_anc_names = ['slope',
                        'aspect']

        config_dict["anc_files"] = {}
        if ascii_loader_bool: # do not sort, use ascii order
            anc_files = setting_dict["input_datasets"]["anc_files"]
        else:
            anc_files = sorted(setting_dict["input_datasets"]["anc_files"])

        if 'topo' in corr_list or 'brdf' in corr_list:
            if not len(images)==len(anc_files):
                feedback_message+=['// Reflectance images and ancillary images do not match in number.']
            elif len(anc_files)==0:
                feedback_message+=['// Ancillary images required for TOPO or BRDF correction, but none is provided.']
            else:
                for i,image in enumerate(images):
                    config_dict["anc_files"][image] = dict(zip(external_anc_names,
                                                          [[anc_files[i],a] for a in range(len(external_anc_names))]))


    config_dict['output_dir'] = setting_dict["export_settings"]["target_export_directory"]+"/"
    config_dict['export_type'] = setting_dict["processing_pipeline"]["export_type"]
    config_dict["use_glt"] = setting_dict["input_datasets"]["use_glt"]

    config_dict["glt_files"] = {}
    if config_dict["use_glt"]:

        if ascii_loader_bool: # do not sort, use ascii order
            glt_files = setting_dict["input_datasets"]["glt_files"]
        else:
            glt_files =  sorted(setting_dict["input_datasets"]["glt_files"])

        if not len(images)==len(glt_files):
            feedback_message+=['// Reflectance images and GLT files do not match.']
        elif len(glt_files)==0:
            feedback_message+=['// GLT files required, but none is provided.']
        else:
            for i,image in enumerate(images):
                config_dict['glt_files'][image] = dict(zip(["glt_x","glt_y"],
                                                [[glt_files[i],a] for a in range(2)]))

    config_dict["corrections"] = corr_list
    config_dict["chunk_size"] = setting_dict["processing_pipeline"]["chunk_size"]

    config_dict["topo"] =  {}
    if ascii_loader_bool: # do not sort, use ascii order
        topo_json_files = setting_dict["input_datasets"]["precomputed"]["topo_json_files"]
    else:
        topo_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["topo_json_files"])

    if 'topo' in corr_list:
        if not len(images)==len(topo_json_files):
            feedback_message+=['// Reflectance images and TOPO coefficient files do not match.']
        elif len(topo_json_files)==0:
            feedback_message+=['// Precomputed TOPO coefficient files required, but none is provided.']
        else:
            config_dict["topo"] = dict(zip(images,topo_json_files))


    config_dict["brdf"]  = {}
    if ascii_loader_bool: # do not sort, use ascii order
        brdf_json_files = setting_dict["input_datasets"]["precomputed"]["brdf_json_files"]
    else:
        brdf_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["brdf_json_files"])

    if 'brdf' in corr_list:
        if not len(images)==len(brdf_json_files):
            feedback_message+=['// Reflectance images and BRDF coefficient files do not match.']
        elif len(brdf_json_files)==0:
            feedback_message+=['// Precomputed BRDF coefficient files required, but none is provided.']
        else:
            config_dict["brdf"] = dict(zip(images,brdf_json_files))

    config_dict["resampling"]  = {}
    if not setting_dict["processing_pipeline"]["resample_method"]=='None':        
        config_dict["resampling"]['type'] = setting_dict["processing_pipeline"]["resample_method"] #'linear' # 'cubic' 'nearest'

    config_dict["masks"] = [["ndi", {'band_1': 850,'band_2': 660,
                                  'min': 0.1,'max': 1.0}]]

    if config_dict["file_type"] =="neon":
        config_dict["masks"]+=[['neon_edge',{'radius': 30}]]                        

    config_dict["trait_models"] = setting_dict["input_datasets"]["trait_files"]

    return config_dict, feedback_message

def process_and_reorder_config(raw_data):
    """
    Internal Python transformation function.
    Cleans, standardizes, and reorders input parameters before JSON serialization.
    """
    inputs = raw_data.get('inputs', {})
    params = raw_data.get('parameters', {})

    bad_band_range = {
        "badband1":[[300,400],[1337,1430],[1800,1960],[2450,2600]],
        "true":[[300,400],[750,2600]],
        "ir":[[300,750],[1337,1430],[1800,1960],[2450,2600]]
    }

    # Standardized & Reordered Schema
    structured_config_from_gui = {
        "job_metadata": {
            "config_version": "1.0.0",
            "created_at": raw_data.get('timestamp', datetime.now(timezone.utc).isoformat() + "Z"),

        },
        "export_settings": {
            "target_export_directory": raw_data.get('export_directory', ''),
        },
        "input_datasets": {

            "refl_images": inputs.get('images', []),
            "total_primary_files": len(inputs.get('images', [])),

            "ascii_loader_enabled": inputs.get('enable_ascii_loader', False),
            "anc_files": inputs.get('anc_files', []),
            "trait_files": inputs.get('trait_files', []),
            "use_glt": inputs.get('use_glt', False),
            "glt_files": inputs.get('glt_files', []),
            "precomputed_enabled": inputs.get('enable_precomputed_correction', False),
            "precomputed": {
                "topo_json_files": inputs.get('topo_json_files', []),
                "brdf_json_files": inputs.get('brdf_json_files', [])
            }
        },
        "processing_pipeline": {
            "data_type": params.get('data_type'),
            "export_type": params.get('data_type_export'),
            "corrections": params.get('corrections', []),
            "chunk_size": params.get('chunk_size',[]),
            "bad_bands": bad_band_range[params.get('bad_bands')],
            "resample_method": params.get('resample_method'),
        }
    }
    structured_config, feedback_message = generate_config(structured_config_from_gui)
    return structured_config, feedback_message


class ScrollableFrame(ttk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_window = ttk.Frame(canvas)

        self.scrollable_window.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=self.scrollable_window, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Mouse wheel scrolling
        self.scrollable_window.bind('<Enter>', self._bound_to_mousewheel)
        self.scrollable_window.bind('<Leave>', self._unbound_to_mousewheel)
        self.canvas = canvas

    def _bound_to_mousewheel(self, event):
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)

    def _unbound_to_mousewheel(self, event):
        self.canvas.unbind_all("<MouseWheel>")

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Local Server Job Configurator for Trait Mapping by HyTools")
        self.geometry("950x800")
        self.configure(bg="#f4f4f5")

        style = ttk.Style(self)
        style.theme_use('clam')

        # Main Scrollable Container
        self.main_frame = ScrollableFrame(self)
        self.main_frame.pack(fill="both", expand=True, padx=10, pady=10)
        self.container = self.main_frame.scrollable_window

        # Variables
        self.primary_files = []
        self.anc_files = []
        self.trait_files = []
        self.topo_files = []
        self.brdf_files = []
        self.glt_files = []

        self.parsed_ascii_data = None

        self.build_ui()

    def create_scrollable_listbox(self, parent, height, pady_val):
        """Helper to create a Listbox with vertical & horizontal scrollbars"""
        container = ttk.Frame(parent)
        container.pack(fill="x", pady=pady_val)

        lb = tk.Listbox(container, height=height, bg="#1f2937", fg="#a7f3d0", selectbackground="#374151")

        y_scroll = ttk.Scrollbar(container, orient="vertical", command=lb.yview)
        x_scroll = ttk.Scrollbar(container, orient="horizontal", command=lb.xview)

        lb.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)

        lb.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")

        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)

        return lb

    def build_ui(self):
        # 1. Primary Files
        f1 = ttk.LabelFrame(self.container, text=" 1. Input Reflectance Images (Required) ", padding=15)
        f1.pack(fill="x", pady=5, padx=5)
        btn_frame_ref = ttk.Frame(f1)
        btn_frame_ref.pack(fill="x", pady=2)
        ttk.Button(btn_frame_ref, text="Browse Reflectance Files...", command=lambda: self.select_files(self.lb_primary, self.primary_files, update_groups=True)).pack(side="left", padx=2)
        ttk.Button(btn_frame_ref, text="Clear Files", command=self.clear_primary_files).pack(side="left", padx=2)
        self.lb_primary = self.create_scrollable_listbox(f1, height=4, pady_val=(0, 10))


        ####
        # ASCII Auto-Loader Frame
        ascii_frame = ttk.LabelFrame(f1, text=" ASCII Spreadsheet Auto-Loader (*.csv, *.txt, *.tsv) ", padding=8)
        ascii_frame.pack(fill="x", pady=5)

        self.ascii_enabled_var = tk.BooleanVar(value=False)
        cb = ttk.Checkbutton(ascii_frame, text="Enable ASCII Spreadsheet Auto-Loader", 
                             variable=self.ascii_enabled_var, command=self.toggle_ascii_ui)
        cb.pack(anchor="w")

        self.ascii_controls_frame = ttk.Frame(ascii_frame)
        self.ascii_controls_frame.pack(fill="x", pady=5)

        ttk.Button(self.ascii_controls_frame, text="Browse & Parse ASCII File...", 
                   command=self.load_ascii_file).pack(anchor="w", pady=2)

        ttk.Label(self.ascii_controls_frame, text="File Preview (First 3 Lines):", font=('TkDefaultFont', 9, 'bold')).pack(anchor="w", pady=(5,0))
        self.ascii_preview_text = tk.Text(self.ascii_controls_frame, height=3, state="disabled", bg="#1e293b", fg="#f1f5f9")
        self.ascii_preview_text.pack(fill="x", pady=2)

        cols_frame = ttk.Frame(self.ascii_controls_frame)
        cols_frame.pack(fill="x", pady=5)

        ttk.Label(cols_frame, text="Reflectance Column:").grid(row=0, column=0, sticky="w", padx=4)
        self.ascii_col_ref = ttk.Combobox(cols_frame, state="readonly")
        self.ascii_col_ref.grid(row=1, column=0, sticky="ew", padx=4)

        ttk.Label(cols_frame, text="Ancillary Column:").grid(row=0, column=1, sticky="w", padx=4)
        self.ascii_col_anc = ttk.Combobox(cols_frame, state="readonly")
        self.ascii_col_anc.grid(row=1, column=1, sticky="ew", padx=4)

        ttk.Label(cols_frame, text="GLT Column:").grid(row=0, column=2, sticky="w", padx=4)
        self.ascii_col_glt = ttk.Combobox(cols_frame, state="readonly")
        self.ascii_col_glt.grid(row=1, column=2, sticky="ew", padx=4)

        ttk.Label(cols_frame, text="TOPO coeff Column:").grid(row=2, column=0, sticky="w", padx=4)
        self.ascii_col_topo = ttk.Combobox(cols_frame, state="readonly")
        self.ascii_col_topo.grid(row=3, column=0, sticky="ew", padx=4)

        ttk.Label(cols_frame, text="BRDF coeff Column:").grid(row=2, column=1, sticky="w", padx=4)
        self.ascii_col_brdf = ttk.Combobox(cols_frame, state="readonly")
        self.ascii_col_brdf.grid(row=3, column=1, sticky="ew", padx=4)


        cols_frame.columnconfigure(0, weight=1)
        cols_frame.columnconfigure(1, weight=1)
        cols_frame.columnconfigure(2, weight=1)
        cols_frame.columnconfigure(3, weight=1)
        cols_frame.columnconfigure(4, weight=1)

        self.ascii_col_ref.bind("<<ComboboxSelected>>", lambda e: self.apply_ascii_data())
        self.ascii_col_anc.bind("<<ComboboxSelected>>", lambda e: self.apply_ascii_data())
        self.ascii_col_glt.bind("<<ComboboxSelected>>", lambda e: self.apply_ascii_data())
        self.ascii_col_topo.bind("<<ComboboxSelected>>", lambda e: self.apply_ascii_data())
        self.ascii_col_brdf.bind("<<ComboboxSelected>>", lambda e: self.apply_ascii_data())

        self.ascii_status_lbl = ttk.Label(self.ascii_controls_frame, text="", foreground="green", font=('TkDefaultFont', 9, 'bold'))
        self.ascii_status_lbl.pack(anchor="w", pady=2)

        self.toggle_ascii_ui()
        ####

        # 2. Ancillary Data
        f3 = ttk.LabelFrame(self.container, text=" 2. (Optional, if no correction is needed) Ancillary images in ENVI format paired with reflectance images (e.g. *obs_ort, *obs ) \n     Leave it blank for NEON AOP images", padding=15)
        f3.pack(fill="x", pady=5, padx=5)
        self.lb_ref = self.create_scrollable_listbox(f3, height=3, pady_val=(0, 10))
        ttk.Button(f3, text="Browse Ancillary Files...", command=lambda: self.select_files(self.lb_ref, self.anc_files)).pack(side="left", padx=2)
        ttk.Button(f3, text="Clear Files", command=lambda: self.clear_listbox(self.lb_ref, self.anc_files)).pack(side="left", padx=2)

        # 3. Precomputed correction model coefficients 
        f4 = ttk.LabelFrame(self.container, text=" 3. Precomputed correction model coefficients ", padding=15)
        f4.pack(fill="x", pady=5, padx=5)

        self.supp_enabled_var = tk.BooleanVar(value=False)
        chk_supp = ttk.Checkbutton(f4, text="Use precomputed correction model coefficients", variable=self.supp_enabled_var, command=self.toggle_supplementary)
        chk_supp.pack(anchor="w", pady=(0, 10))

        self.supp_container = ttk.Frame(f4)
        self.supp_container.pack(fill="x")

        ttk.Label(self.supp_container, text="Selector TOPO coefficients (JSON Files paired with reflectance images)").pack(anchor="w")        
        btn_f1 = ttk.Frame(self.supp_container)
        btn_f1.pack(fill="x")
        self.btn_supp_a = ttk.Button(btn_f1, text="Browse TOPO JSON Files...", command=lambda: self.select_files(self.lb_supp_a, self.topo_files, filetypes=[("JSON Files","*.json;*.JSON")]))
        self.btn_supp_a.pack(side="left", pady=(0, 10))
        ttk.Button(btn_f1, text="Clear", command=lambda: self.clear_listbox(self.lb_supp_a, self.topo_files)).pack(side="left", pady=(0, 10),padx=2)
        self.lb_supp_a = self.create_scrollable_listbox(self.supp_container, height=3, pady_val=(0, 5))

        ttk.Label(self.supp_container, text="Selector BRDF coefficients (JSON Files paired with reflectance images)").pack(anchor="w")
        btn_f2 = ttk.Frame(self.supp_container)
        btn_f2.pack(fill="x")        
        self.btn_supp_b = ttk.Button(btn_f2, text="Browse BRDF JSON Files...", command=lambda: self.select_files(self.lb_supp_b, self.brdf_files, filetypes=[("JSON Files","*.json;*.JSON")]))
        self.btn_supp_b.pack(side="left", pady=(0, 10))
        ttk.Button(btn_f2, text="Clear", command=lambda: self.clear_listbox(self.lb_supp_b, self.brdf_files)).pack(side="left", pady=(0, 10),padx=2)
        self.lb_supp_b = self.create_scrollable_listbox(self.supp_container, height=3, pady_val=(0, 5))

        self.toggle_supplementary() # Initialize disabled state

        # 4. Use GLT
        f8 = ttk.LabelFrame(self.container, text=" 4. Use External GLT for location adjustment ", padding=15)
        f8.pack(fill="x", pady=5, padx=5)
        
        self.enabled_glt = tk.BooleanVar(value=False)
        chk_glt = ttk.Checkbutton(f8, text="Use GLT file", variable=self.enabled_glt, command=self.toggle_glt)
        chk_glt.pack(anchor="w", pady=(0, 10))

        self.glt_container = ttk.Frame(f8)
        self.glt_container.pack(fill="x")

        ttk.Label(self.glt_container, text="Selector GLT (2-band ENVI Files paired with reflectance images)").pack(anchor="w")
        self.lb_glt = self.create_scrollable_listbox(self.glt_container, height=3, pady_val=(0, 5))
        self.btn_glt = ttk.Button(self.glt_container, text="Browse GLT Files...", command=lambda: self.select_files(self.lb_glt, self.glt_files))
        self.btn_glt.pack(anchor="w", pady=(0, 10))

        self.toggle_glt() # Initialize disabled state


        # 5. Parameters
        f5 = ttk.LabelFrame(self.container, text=" 5. Correction Settings ", padding=15)
        f5.pack(fill="x", pady=5, padx=5)

        p_grid = ttk.Frame(f5)
        p_grid.pack(fill="x")

        col1 = ttk.Frame(p_grid)
        col1.pack(side="left", fill="both", expand=True,padx=5)
        col2 = ttk.Frame(p_grid)
        col2.pack(side="right", fill="both", expand=True,padx=5)

        self.var_datatype = tk.StringVar(value="envi")
        ttk.Label(col1, text="Data Type for Input Reflectance Images", font=("", 9, "bold")).pack(anchor="w")
        ttk.Radiobutton(col1, text="ENVI", variable=self.var_datatype, value="envi").pack(anchor="w")
        ttk.Radiobutton(col1, text="NEON AOP HDF5 (*.h5)", variable=self.var_datatype, value="neon").pack(anchor="w")
        ttk.Radiobutton(col1, text="AVIRIS NetCDF (*.nc)", variable=self.var_datatype, value="ncav").pack(anchor="w")
        ttk.Radiobutton(col1, text="EMIT NetCDF (*.nc)", variable=self.var_datatype, value="emit").pack(anchor="w")
        ttk.Radiobutton(col1, text="Tanager HDF5 (*.h5)", variable=self.var_datatype, value="tanager").pack(anchor="w", pady=(0, 10))

        self.var_datatype_export = tk.StringVar(value="envi")
        ttk.Label(col1, text="Data Type for Output Reflectance Images", font=("", 9, "bold")).pack(anchor="w")
        ttk.Radiobutton(col1, text="ENVI", variable=self.var_datatype_export, value="envi").pack(anchor="w")
        ttk.Radiobutton(col1, text="NetCDF", variable=self.var_datatype_export, value="netcdf").pack(anchor="w", pady=(0, 10))

        self.var_in_band_subset = tk.StringVar(value="badband1")
        ttk.Label(col1, text="Good band range in reflectance images", font=("", 9, "bold")).pack(anchor="w")
        ttk.Radiobutton(col1, text="Bands outside ([[300,400],[1337,1430],[1800,1960],[2450,2600]])", variable=self.var_in_band_subset, value="badband1").pack(anchor="w")
        ttk.Radiobutton(col1, text="Visible bands (400-750)", variable=self.var_in_band_subset, value="true").pack(anchor="w")
        ttk.Radiobutton(col1, text="Bands outside ([[300,750],[1337,1430],[1800,1960],[2450,2600]])", variable=self.var_in_band_subset, value="ir").pack(anchor="w", pady=(0, 10))


        self.var_topo_bool= tk.StringVar(value="")
        self.var_brdf_bool = tk.StringVar(value="")
        self.var_glint_bool = tk.StringVar(value="")        
        ttk.Label(col2, text="Correction", font=("", 9, "bold")).pack(anchor="w")
        ttk.Checkbutton(col2, text="Topographic Correction (topo)", variable=self.var_topo_bool, onvalue="topo", offvalue="").pack(anchor="w") #,command=self.update_suffix_entry
        ttk.Checkbutton(col2, text="BRDF correction (brdf)", variable=self.var_brdf_bool, onvalue="brdf", offvalue="").pack(anchor="w") # ,command=self.update_suffix_entry


        chunk_col_options = {
                "Full Column":"col",
                "16":16,
                "64":64,
                "128":128,
                "256":256,
        }
        chunk_row_options = {
                "Full Row":"row",
                "16":16,
                "64":64,
                "128":128,
                "256":256,
        }
        self.chunk_row_size = 64
        self.chunk_col_size = 64

        # Create the Combobox
        ttk.Label(col2, text="Select Chunk Row Size", font=("", 9, "bold")).pack(anchor="w")
        self.row_combo = ttk.Combobox(col2, values=list(chunk_row_options.keys()), state="readonly")
        self.row_combo.set("Select row size")  # Set placeholder text
        self.row_combo.pack(anchor="w", pady=(0, 10))
        self.row_combo.current(2)


        ttk.Label(col2, text="Select Chunk Column Size", font=("", 9, "bold")).pack(anchor="w")
        self.col_combo = ttk.Combobox(col2, values=list(chunk_col_options.keys()), state="readonly")
        self.col_combo.set("Select column size")  # Set placeholder text
        self.col_combo.pack(anchor="w", pady=(0, 10))
        self.col_combo.current(2)

        self.chunk_combination = ttk.Label(col2, text=f"Chunk Size [{chunk_row_options[self.row_combo.get()]} × {chunk_col_options[self.col_combo.get()]}]", font=("", 9))  #{chunk_row_options[self.row_combo.get()],chunk_col_options[self.col_combo]}
        self.chunk_combination.pack(anchor="w")

        self.row_combo.bind("<<ComboboxSelected>>", lambda event: self.on_chunk_selected(event, chunk_row_options))
        self.col_combo.bind("<<ComboboxSelected>>",  lambda event: self.on_chunk_selected(event, chunk_col_options))

        resample_list = ['None','nearest','linear','cubic']
        ttk.Label(col2, text="Select Resampling Method", font=("", 9, "bold")).pack(anchor="w")
        self.col_resample = ttk.Combobox(col2, values=resample_list, state="readonly")
        self.col_resample.set("Select Resampling Method")  # Set placeholder text
        self.col_resample.pack(anchor="w", pady=(0, 10))
        self.col_resample.current(1)

        # Bind the selection event to the trigger function

        # 6. Trait file loading
        f2 = ttk.LabelFrame(self.container, text=" 6. Select Trait Model JSON File ", padding=15)
        f2.pack(fill="x", pady=5, padx=5)

        self.lb_trait = self.create_scrollable_listbox(f2, height=3, pady_val=(0, 5))
        self.btn_trait = ttk.Button(f2, text="Browse Trait JSON Files...", command=lambda: self.select_files(self.lb_trait, self.trait_files, filetypes=[("JSON Files","*.json;*.JSON")]))
        self.btn_trait.pack(side="left", pady=(0, 10))
        ttk.Button(f2, text="Clear", command=lambda: self.clear_listbox(self.lb_trait, self.trait_files)).pack(side="left", pady=(0, 10),padx=2)

        # 7. Export Settings
        f6 = ttk.LabelFrame(self.container, text=" 7. Export Folder Setup (Stored inside JSON Parameters) ", padding=15)
        f6.pack(fill="x", pady=5, padx=5)

        export_grid = ttk.Frame(f6)
        export_grid.pack(fill="x", pady=5, padx=5)

        self.export_dir_var = tk.StringVar(value="Not set")
        ttk.Label(export_grid, textvariable=self.export_dir_var, foreground="#d946ef", font=("Courier", 10)).grid(row=0, column=0, sticky="w", pady=2,padx=5)
        ttk.Button(export_grid, text="Select Export Directory...", command=lambda: self.select_directory(self.export_dir_var)).grid(row=1, column=0, sticky="w", pady=2,padx=5)

        # 8. Save Settings
        f7 = ttk.LabelFrame(self.container, text=" 8. Save JSON Configuration File ", padding=15)
        f7.pack(fill="x", pady=5, padx=5)
        
        save_grid = ttk.Frame(f7)
        save_grid.pack(fill="x")

        ttk.Label(save_grid, text="JSON Save Folder:", font=("", 9, "bold")).grid(row=0, column=0, sticky="w", pady=2)
        self.json_save_dir_var = tk.StringVar(value=os.path.abspath(os.path.dirname(__file__)))
        ttk.Label(save_grid, textvariable=self.json_save_dir_var, foreground="#d946ef", font=("Courier", 10)).grid(row=1, column=0, sticky="w", pady=(0, 10), padx=(0, 20))
        ttk.Button(save_grid, text="Select JSON Save Folder...", command=lambda: self.select_directory(self.json_save_dir_var)).grid(row=2, column=0, sticky="w")

        ttk.Label(save_grid, text="Configuration File Name:", font=("", 9, "bold")).grid(row=0, column=1, sticky="w", pady=2)
        self.json_filename_var = tk.StringVar(value="local_job_config.json")
        ttk.Entry(save_grid, textvariable=self.json_filename_var, width=30).grid(row=1, column=1, sticky="nw")

        # Action Buttons & Preview
        btn_frame = ttk.Frame(self.container)
        btn_frame.pack(fill="x", pady=15, padx=5)
        ttk.Button(btn_frame, text="1. Generate & Preview JSON", command=self.generate_preview).pack(side="left", expand=True, fill="x", padx=(0,5))
        ttk.Button(btn_frame, text="2. Save JSON to Disk", command=self.save_json).pack(side="right", expand=True, fill="x", padx=(5,0))

        # JSON Preview Area
        preview_frame = ttk.Frame(self.container)
        preview_frame.pack(fill="both", expand=True, pady=5, padx=5)

        self.txt_preview = tk.Text(preview_frame, height=12, bg="#111827", fg="#38bdf8", font=("Courier", 10), wrap="none")

        y_scroll = ttk.Scrollbar(preview_frame, orient="vertical", command=self.txt_preview.yview)
        x_scroll = ttk.Scrollbar(preview_frame, orient="horizontal", command=self.txt_preview.xview)

        self.txt_preview.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)

        self.txt_preview.grid(row=0, column=0, sticky="nsew")
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll.grid(row=1, column=0, sticky="ew")

        preview_frame.rowconfigure(0, weight=1)
        preview_frame.columnconfigure(0, weight=1)

    def select_files(self, listbox, target_list, update_groups=False, filetypes=None):
        if filetypes is None:
            files = filedialog.askopenfilenames(title="Select Files")
        else:
            files = filedialog.askopenfilenames(title="Select Files", filetypes=filetypes)
        if files:
            target_list.clear()
            listbox.delete(0, tk.END)
            for f in files:
                target_list.append(f)
                listbox.insert(tk.END, f)

    def refresh_listbox(self, listbox, items):
        listbox.delete(0, tk.END)
        for item in items:
            listbox.insert(tk.END, item)

    def select_directory(self, string_var):
        dir_path = filedialog.askdirectory(title="Select Directory")
        if dir_path:
            string_var.set(dir_path)

    ###ascii
    def toggle_ascii_ui(self):
        state = "normal" if self.ascii_enabled_var.get() else "disabled"
        for child in self.ascii_controls_frame.winfo_children():
            if isinstance(child, (ttk.Frame, ttk.LabelFrame)):
                for sub in child.winfo_children():
                    sub.configure(state=state)
            else:
                child.configure(state=state)

    # ------------------ ASCII Loader Logic ------------------
    def load_ascii_file(self):
        filepath = filedialog.askopenfilename(
            title="Select ASCII Spreadsheet File",
            filetypes=[("Spreadsheet files", "*.csv *.txt *.tsv *.dat"), ("All files", "*.*")]
        )
        if not filepath:
            return

        try:
            with open(filepath, 'r', encoding='utf-8', errors='strict') as f:
                lines = [f.readline() for _ in range(50)]

            raw_preview = [l.strip('\r\n') for l in lines[:3] if l.strip()]
            if not raw_preview:
                messagebox.showerror("Error", "Selected file is empty.")
                return

            header_line = raw_preview[0]
            if ',' in header_line: delimiter = ','
            elif '\t' in header_line: delimiter = '\t'
            elif ';' in header_line: delimiter = ';'
            else: delimiter = None

            with open(filepath, 'r', encoding='utf-8') as f:
                if delimiter:
                    reader = csv.reader(f, delimiter=delimiter)
                else:
                    reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
                all_rows = [row for row in reader if row and not row[0].startswith('#')]

            if not all_rows:
                messagebox.showerror("Error", "No valid data rows found in ASCII file.")
                return

            headers = [h.strip() for h in all_rows[0]]
            data_rows = [[cell.strip() for cell in row] for row in all_rows[1:]]

            self.parsed_ascii_data = {'headers': headers, 'rows': data_rows}

            # Render Preview Text
            self.ascii_preview_text.config(state="normal")
            self.ascii_preview_text.delete('1.0', tk.END)
            self.ascii_preview_text.insert(tk.END, "\n".join(raw_preview))
            self.ascii_preview_text.config(state="disabled")

            # Update Dropdown Columns
            cols = [f"{hdr} (Col {idx+1})" for idx, hdr in enumerate(headers)]
            self.ascii_col_ref['values'] = cols
            self.ascii_col_anc['values'] = cols
            self.ascii_col_glt['values'] = cols
            self.ascii_col_topo['values'] = cols
            self.ascii_col_brdf['values'] = cols

            # Intelligent Default Mapping
            for idx, hdr in enumerate(headers):
                lower = hdr.lower()
                if any(k in lower for k in ['ref', 'primary', 'image', 'reflectance']):
                    self.ascii_col_ref.current(idx)
                if any(k in lower for k in ['anc', 'aux', 'obs', 'obs_ort']):
                    self.ascii_col_anc.current(idx)
                if any(k in lower for k in ['glt', 'lookup']):
                    self.ascii_col_glt.current(idx)
                if any(k in lower for k in ['topo']):
                    self.ascii_col_topo.current(idx)
                if any(k in lower for k in ['brdf']):
                    self.ascii_col_brdf.current(idx)

            if self.ascii_col_ref.current() < 0: self.ascii_col_ref.current(0)
            if self.ascii_col_anc.current() < 0: self.ascii_col_anc.current(0)
            if self.ascii_col_glt.current() < 0: self.ascii_col_glt.current(0)
            if self.ascii_col_topo.current() < 0: self.ascii_col_topo.current(0)
            if self.ascii_col_brdf.current() < 0: self.ascii_col_brdf.current(0)

            self.apply_ascii_data()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse ASCII file: {str(e)}")

    def apply_ascii_data(self):
        if not self.parsed_ascii_data or not self.parsed_ascii_data['rows']:
            return

        ref_idx = self.ascii_col_ref.current()
        anc_idx = self.ascii_col_anc.current()
        glt_idx = self.ascii_col_glt.current()
        brdf_idx = self.ascii_col_brdf.current()
        topo_idx = self.ascii_col_topo.current()

        if ref_idx < 0 or anc_idx < 0 or glt_idx<0 or brdf_idx<0 or topo_idx < 0:
            return

        rows = self.parsed_ascii_data['rows']

        self.primary_files.clear()
        self.anc_files.clear()
        self.glt_files.clear()
        self.topo_files.clear()
        self.brdf_files.clear()

        for row in rows:
            if ref_idx < len(row) and row[ref_idx]:
                ref_file = row[ref_idx]
                #if ref_file not in self.primary_files:
                #    self.primary_files.append(ref_file)
                self.primary_files.append(ref_file)
                #subgroup = row[topo_idx] if (topo_idx < len(row) and row[topo_idx]) else 'Unassigned'
                #if subgroup not in unique_subgroups:
                #    unique_subgroups.append(subgroup)
                #file_to_subgroup[ref_file] = subgroup

            if anc_idx < len(row) and row[anc_idx]:
                ancil_file = row[anc_idx]
                #if ancil_file not in self.anc_files:
                #    self.anc_files.append(ancil_file)
                self.anc_files.append(ancil_file)

            if glt_idx < len(row) and row[glt_idx]:
                glt_file = row[glt_idx]
                #if glt_file not in self.glt_files:
                #    self.glt_files.append(glt_file)
                self.glt_files.append(glt_file)

            if topo_idx < len(row) and row[topo_idx]:
                topo_file = row[topo_idx]
                #if topo_file not in self.topo_files:
                #    self.topo_files.append(topo_file)
                self.topo_files.append(topo_file)

            if brdf_idx < len(row) and row[brdf_idx]:
                brdf_file = row[brdf_idx]
                #if brdf_file not in self.brdf_files:
                #    self.brdf_files.append(brdf_file)
                self.brdf_files.append(brdf_file)
    
        self.refresh_listbox(self.lb_primary, self.primary_files)
        self.refresh_listbox(self.lb_ref, self.anc_files)
        self.refresh_listbox(self.lb_glt, self.glt_files)
        self.refresh_listbox(self.lb_supp_a, self.topo_files)
        self.refresh_listbox(self.lb_supp_b, self.brdf_files)

        self.ascii_status_lbl.config(
            text=f"✓ Auto-loaded {len(self.primary_files)} Reflectance, {len(self.anc_files)} Ancillary files."
        )


    def toggle_supplementary(self):
        state = 'normal' if self.supp_enabled_var.get() else 'disabled'
        self.lb_supp_a.configure(state=state)
        self.btn_supp_a.configure(state=state)
        self.lb_supp_b.configure(state=state)
        self.btn_supp_b.configure(state=state)

    def toggle_glt(self):
        state = 'normal' if self.enabled_glt.get() else 'disabled'
        self.lb_glt.configure(state=state)
        self.btn_glt.configure(state=state)

    def update_grouping_ui(self, *args):
        try:
            num = int(self.num_groups_var.get())
            if num < 1: num = 1
        except ValueError:
            num = 2

        # Adjust tracking lists
        while len(self.group_name_vars) < num:
            self.group_name_vars.append(tk.StringVar(value=f"Group {len(self.group_name_vars)+1}"))
        if len(self.group_name_vars) > num:
            self.group_name_vars = self.group_name_vars[:num]

        # Redraw Group Names

        for widget in self.group_names_frame.winfo_children():
            widget.destroy()
        print("self.group_name_vars",self.group_name_vars)
        for i, var in enumerate(self.group_name_vars):
            row = ttk.Frame(self.group_names_frame)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"{i+1}.", width=3).pack(side="left")
            entry = ttk.Entry(row, textvariable=var)
            entry.pack(side="left", fill="x", expand=True)
            var.trace_add("write", lambda *args: self.refresh_combobox_options())

        self.update_file_assignments()

    ######

    def clear_primary_files(self):
        self.primary_files.clear()
        self.lb_primary.delete(0, tk.END)


    def clear_listbox(self, listbox, target_list):
        target_list.clear()
        listbox.delete(0, tk.END)

    def on_chunk_selected(self, event,  lookup_dict):
        """Callback function triggered when a user selects an item."""
        # Get the selected key (Label) from the combobox
        selected_label = event.widget.get()
        #widget_id = event.widget

        # Fetch the corresponding backend ID (Value) from the dictionary
        associated_value = lookup_dict[selected_label]

        # Update the display labels with the results

        if event.widget==self.col_combo:
            self.chunk_col_size = associated_value
        elif event.widget==self.row_combo:
            self.chunk_row_size = associated_value

        self.chunk_combination.config(text=f"Chunk Size = [{self.chunk_row_size},{self.chunk_col_size}]")


    def collect_raw_payload(self):

        preprocs = [v for v in [self.var_topo_bool.get(), self.var_brdf_bool.get()] if v]

        return {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "export_directory": self.export_dir_var.get() if self.export_dir_var.get() != "Not set" else "",
            "inputs": {
                "images":self.primary_files,

                "enable_ascii_loader": self.ascii_enabled_var.get(),
                "anc_files": self.anc_files,
                "trait_files": self.trait_files,
                "enable_precomputed_correction": self.supp_enabled_var.get(),
                "topo_json_files": self.topo_files if self.supp_enabled_var.get() else [], #topo
                "brdf_json_files": self.brdf_files if self.supp_enabled_var.get() else [], #brdf
                "use_glt": self.enabled_glt.get(),
                "glt_files": self.glt_files if self.enabled_glt.get() else [],                
            },
            "parameters": {
                "data_type": self.var_datatype.get(),
                "data_type_export": self.var_datatype_export.get(),
                "bad_bands":self.var_in_band_subset.get(),

                "corrections": preprocs,
                "chunk_size": [self.chunk_row_size, self.chunk_col_size],
                "resample_method": self.col_resample.get(),
            }
        }

    def generate_preview(self):
        raw_payload = self.collect_raw_payload()
        processed_config, feedback_message = process_and_reorder_config(raw_payload)

        self.txt_preview.delete("1.0", tk.END)
        self.txt_preview.insert(tk.END, "// Processed & Reordered JSON Config:\n\n")
        if len(feedback_message)>0:
            for msg in feedback_message:
                self.txt_preview.insert(tk.END, msg+"\n\n")
        self.txt_preview.insert(tk.END, json.dumps(processed_config, indent=4))
        self.main_frame.canvas.yview_moveto(1.0)

    def save_json(self):
        raw_payload = self.collect_raw_payload()
        processed_config, feedback_message = process_and_reorder_config(raw_payload)

        if len(feedback_message)>0:
            for msg in feedback_message:
                self.txt_preview.insert(tk.END, msg+"\n\n")
            self.txt_preview.delete("1.0", tk.END)
            self.txt_preview.insert(tk.END, f"// File is not saved due to warnings.\n\n")
            messagebox.showerror("Error", f"Failed to save file:\n")
            return

        save_dir = self.json_save_dir_var.get()
        filename = self.json_filename_var.get().strip()
        if not filename.endswith('.json'):
            filename += '.json'

        full_path = os.path.join(save_dir, filename)

        try:
            with open(full_path, 'w') as f:
                json.dump(processed_config, f, indent=4)

            self.txt_preview.delete("1.0", tk.END)
            self.txt_preview.insert(tk.END, f"// SUCCESS: File written to disk!\n// Destination: {full_path}\n\n")
            self.txt_preview.insert(tk.END, json.dumps(processed_config, indent=4))
            messagebox.showinfo("Success", f"Configuration successfully saved to:\n{full_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file:\n{str(e)}")

if __name__ == "__main__":
    app = App()
    app.mainloop()

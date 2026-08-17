from flask import Flask, request, jsonify, render_template_string
import os
import sys
import json
import csv
from datetime import datetime, timezone
import numpy as np

app = Flask(__name__)

def generate_config_v2(setting_dict):


    feedback_message = []
    corr_list = setting_dict["processing_pipeline"]["corrections"]
    ascii_loader_bool = setting_dict["input_datasets"]["ascii_loader_enabled"]


    config_dict = {}
    config_dict["job_metadata"] = setting_dict["job_metadata"]
    config_dict["job_metadata"]["warning_message"] = {}
    config_dict['bad_bands'] = setting_dict["processing_pipeline"]["bad_bands"]
    config_dict["file_type"] = setting_dict["processing_pipeline"]["data_type"]

    if ascii_loader_bool: # do not sort, use ascii order
        images = setting_dict["input_datasets"]["refl_images"]
    else:
        images = sorted(setting_dict["input_datasets"]["refl_images"])    

    config_dict["input_files"] = images

    if len(images)<1:
        feedback_message+=['// No reflectance images is provided.']

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
    elif config_dict["file_type"] =="ncav":

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

    config_dict['export'] = {}
    config_dict['export']['coeffs'] = 'coeffs' in setting_dict["processing_pipeline"]["tasks_to_execute"]
    config_dict['export']['image'] = 'image' in setting_dict["processing_pipeline"]["tasks_to_execute"]
    config_dict['export']['masks'] = 'masks' in setting_dict["processing_pipeline"]["tasks_to_execute"]
    config_dict['export']['subset_waves']  = setting_dict["processing_pipeline"]["subset_waves"]
    config_dict['export']['output_dir'] = setting_dict["export_settings"]["target_export_directory"]+"/"
    config_dict['export']['image_format'] = setting_dict["processing_pipeline"]["export_type"]
    config_dict['export']["suffix"] = setting_dict["processing_pipeline"]["suffix"]
    config_dict['export']["use_glt"] = setting_dict["input_datasets"]["use_glt"]

    if config_dict['export']["use_glt"]:
        glt_files =  setting_dict["input_datasets"]["glt_files"]

        #if ascii_loader_bool: # do not sort, use ascii order
        #    glt_files = setting_dict["input_datasets"]["glt_files"]
        #else:
        #    glt_files =  sorted(setting_dict["input_datasets"]["glt_files"])        

        config_dict["glt_files"] = {}
        if not len(images)==len(glt_files):
            feedback_message+=['// Reflectance images and GLT files do not match.']
        elif len(glt_files)==0:
            feedback_message+=['// GLT files required, but none is provided.']
        else:
            for i,image in enumerate(images):
                config_dict['glt_files'][image] = dict(zip(["glt_x","glt_y"],
                                                [[glt_files[i],a] for a in range(2)]))

    config_dict["corrections"] = corr_list

    config_dict["topo"] =  {}
    config_dict["topo"]['type'] =  setting_dict["processing_pipeline"]["topo_method"]
    topo_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["topo_json_files"])
    #if ascii_loader_bool: # do not sort, use ascii order
    #    topo_json_files = setting_dict["input_datasets"]["precomputed"]["topo_json_files"]
    #else:
    #    topo_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["topo_json_files"])

    if config_dict["topo"]['type']=='precomputed':
        if not len(images)==len(topo_json_files):
            feedback_message+=['// Reflectance images and TOPO coefficient files do not match.']
        elif len(topo_json_files)==0:
            feedback_message+=['// Precomputed TOPO coefficient files required, but none is provided.']
        else:
            config_dict["topo"]['coeff_files'] = dict(zip(images,topo_json_files))
    else:
        config_dict["topo"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                    'min': 0.1,'max': 1.0}],
                                            ['ancillary',{'name':'slope',
                                                        'min': np.radians(5),'max':'+inf' }],
                                            ['ancillary',{'name':'cosine_i',
                                                        'min': 0.12,'max':'+inf' }],
                                            ['cloud',{'method':'zhai_2018',
                                                    'cloud':True,'shadow':True,
                                                    'T1': 0.01,'t2': 1/10,'t3': 1/4,
                                                    't4': 1/2,'T7': 9,'T8': 9}]]
        if config_dict["file_type"] =="neon":
            config_dict["topo"]['calc_mask']+=[['neon_edge',{'radius': 30}]]
        config_dict["topo"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                    'min': 0.1,'max': 1.0}],
                                            ['ancillary',{'name':'slope',
                                                        'min': np.radians(5),'max':'+inf' }],
                                            ['ancillary',{'name':'cosine_i',
                                                        'min': 0.12,'max':'+inf' }]]
        config_dict["topo"]['c_fit_type'] = 'nnls'

        if len(setting_dict["input_datasets"]["group_images"])>1:
            config_dict["topo"]["subgrouped"] = True
            config_dict["topo"]["sample_perc"] = 0.2
            config_dict["topo"]["subgroup"] = {}
            for group_name in setting_dict["input_datasets"]["group_images"]:
                for image in setting_dict["input_datasets"]["group_images"][group_name]:
                    config_dict["topo"]["subgroup"][image] = group_name


    config_dict["brdf"]  = {}

    # Options are 'line','scene', or a float for a custom solar zn
    # Custom solar zenith angle should be in radians
    config_dict["brdf"]['solar_zn_type'] = setting_dict["processing_pipeline"]["brdf_solar_zn"]

    config_dict["brdf"]['type'] =  setting_dict["processing_pipeline"]["brdf_method"]

    brdf_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["brdf_json_files"])
    #if ascii_loader_bool: # do not sort, use ascii order
    #    brdf_json_files = setting_dict["input_datasets"]["precomputed"]["brdf_json_files"]
    #else:
    #    brdf_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["brdf_json_files"])

    if config_dict["brdf"]['type']=='precomputed':
        if not len(images)==len(brdf_json_files):
            feedback_message+=['// Reflectance images and BRDF coefficient files do not match.']
        elif len(brdf_json_files)==0:
            feedback_message+=['// Precomputed BRDF coefficient files required, but none is provided.']
        else:
            config_dict["brdf"]['coeff_files'] = dict(zip(images,brdf_json_files))
    else:
        if setting_dict["processing_pipeline"]["brdf_group_mode"]=='group':
            config_dict["brdf"]['grouped'] =  True
        else:
            config_dict["brdf"]['grouped'] =  False

        config_dict["brdf"]['geometric'] = 'li_dense_r'
        config_dict["brdf"]['volume'] = 'ross_thick'
        config_dict["brdf"]["b/r"] = 2.5
        config_dict["brdf"]["h/b"] = 2
        config_dict["brdf"]['sample_perc'] = 0.1
        config_dict["brdf"]['interp_kind'] = 'linear'
        config_dict["brdf"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                    'min': 0.1,'max': 1.0}],
                                            ['kernel_finite',{}],
                                            ['ancillary',{'name':'sensor_zn',
                                                        'min':np.radians(2),'max':'inf' }],
                                            ['cloud',{'method':'zhai_2018',
                                                    'cloud':True,'shadow':True,
                                                    'T1': 0.01,'t2': 1/10,'t3': 1/4,
                                                    't4': 1/2,'T7': 9,'T8': 9}]]
        if config_dict["file_type"] =="neon":
            config_dict["brdf"]['calc_mask']+=[['neon_edge',{'radius': 30}]]
        config_dict["brdf"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                    'min': 0.05,'max': 1.0}]]

        # ## Flex dynamic NDVI params
        config_dict["brdf"]['bin_type'] = 'dynamic'
        config_dict["brdf"]['num_bins'] = 18
        config_dict["brdf"]['ndvi_bin_min'] = 0.05
        config_dict["brdf"]['ndvi_bin_max'] = 1.0
        config_dict["brdf"]['ndvi_perc_min'] = 10
        config_dict["brdf"]['ndvi_perc_max'] = 95

    config_dict["glint"]  = {}
    config_dict['glint']['type'] = setting_dict["processing_pipeline"]["glint_method"]
    config_dict['glint']['correction_wave'] = setting_dict["processing_pipeline"]["glint_ref_wave"]
    config_dict["glint"]['apply_mask'] = [["ndi", {'band_1': 550,'band_2': 2150,
                                                'min': 0, 'max': 1}],
                                          ["ndi", {'band_1': 850,'band_2': 660,
                                                'min': -1,'max': 0.1}]]
    config_dict["resample"]  = False
    config_dict['num_cpus'] = len(images)

    config_dict["job_metadata"]["warning_message"] = feedback_message

    return config_dict

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
    export_band_subset = {
        "full":[],
        "truecolor": [660,550,440],
        "subset1": [440,550,660,850,976,1650,2217]
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
            "total_primary_files": sum(len(files) for files in inputs.get('primary_files_grouped', {}).values()),
            "refl_images": inputs.get('images', []),
            "group_images": inputs.get('primary_files_grouped', {}),
            "anc_files": inputs.get('anc_files', []),
            "use_glt": inputs.get('use_glt', False),
            "glt_files": inputs.get('glt_files', []),
            "precomputed_enabled": inputs.get('enable_precomputed_correction', False),
            "precomputed": {
                "topo_json_files": inputs.get('topo_json_files', []),
                "brdf_json_files": inputs.get('brdf_json_files', [])
            },
            "ascii_loader_enabled": inputs.get('enable_ascii_loader', False),
        },
        "processing_pipeline": {
            "data_type": params.get('data_type'),
            "export_type": params.get('data_type_export'),
            "tasks_to_execute": params.get('tasks', []),
            "corrections": params.get('corrections', []),
            "bad_bands":bad_band_range[params.get('bad_bands')],
            "subset_waves":export_band_subset[params.get('subset_waves')],
            "suffix":params.get('suffix'),
            "topo_method":params.get('topo_method'),
            "brdf_method":params.get('brdf_method'),
            "brdf_group_mode":params.get('brdf_group_mode'),
            "brdf_solar_zn":params.get('brdf_solar_zn'),
            "glint_method":params.get('glint_method'),
            "glint_ref_wave":params.get('glint_ref_wave'),
        }
    }
    structured_config = generate_config_v2(structured_config_from_gui)
    return structured_config


HTML_TEMPLATE = r"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Server Job Configurator for Image Correction</title>
  <style>
    body { font-family: system-ui, -apple-system, sans-serif; background: #f4f4f5; padding: 20px; color: #333; max-width: 950px; margin: 0 auto; }
    .card { background: white; padding: 24px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); margin-bottom: 20px; }
    h3 { margin-top: 0; border-bottom: 1px solid #eee; padding-bottom: 8px; color: #111827; }
    h4 { margin-bottom: 8px; color: #374151; font-size: 15px; }
    
    .browser-header { background: #e5e7eb; padding: 10px; border-radius: 4px; font-family: monospace; font-size: 13px; margin-bottom: 10px; word-break: break-all; display: flex; justify-content: space-between; align-items: center; }
    .file-list { max-height: 180px; overflow-y: auto; border: 1px solid #ccc; border-radius: 4px; padding: 0; margin: 0; list-style: none; background: #fff;}
    .file-list li { padding: 6px 12px; border-bottom: 1px solid #f0f0f0; display: flex; align-items: center; }
    .file-list li:hover { background: #f9fafb; }
    .dir-link { color: #2563eb; cursor: pointer; font-weight: 600; text-decoration: none; display: block; width: 100%; user-select: none; }
    
    .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
    .grid-3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 15px; }
    .form-group { margin-bottom: 16px; background: #fafafa; padding: 12px; border-radius: 6px; border: 1px solid #e5e7eb; }
    .options-group { display: flex; flex-direction: column; gap: 8px; margin-top: 8px; }
    
    input[type="text"], input[type="number"], select { padding: 6px 10px; border: 1px solid #ccc; border-radius: 4px; font-size: 14px; }
    
    .btn-group { display: flex; gap: 12px; margin-top: 15px; }
    button { background: #2563eb; color: white; border: none; padding: 12px 20px; border-radius: 6px; cursor: pointer; font-weight: bold; font-size: 15px; flex: 1; }
    button:hover { background: #1d4ed8; }
    .btn-secondary { background: #10b981; }
    .btn-secondary:hover { background: #059669; }
    .btn-save { background: #8b5cf6; }
    .btn-save:hover { background: #7c3aed; }
    
    .selected-box { 
      background: #1f2937; 
      color: #a7f3d0; 
      padding: 10px; 
      border-radius: 4px; 
      font-family: monospace; 
      font-size: 12px; 
      max-height: 110px; 
      overflow-y: auto; 
      margin-bottom: 10px;
      border: 1px solid #374151;
    }
    
    .path-highlight { color: #d946ef; font-weight: bold; font-family: monospace; }
    
    pre#json-preview { 
      background: #111827; 
      color: #38bdf8; 
      padding: 16px; 
      border-radius: 6px; 
      overflow-x: auto; 
      overflow-y: auto; 
      max-height: 250px; 
      line-height: 1.4;
      font-size: 13px;
      display: none; 
      margin-top: 15px; 
      border: 1px solid #374151;
    }
    
    .grouping-table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px;}
    .grouping-table th, .grouping-table td { padding: 8px; border: 1px solid #e5e7eb; text-align: left; }
    .grouping-table th { background: #f9fafb; font-weight: 600; }

    .disabled-zone {
      pointer-events: none;
      opacity: 0.4;
      user-select: none;
      filter: grayscale(1);
    }
  </style>
</head>
<body>

  <!-- FILE SELECTOR 1 -->
  <div class="card">
    <h3>1. Input Reflectance Images (Required)</h3>
    <div class="browser-header"><span id="b1-current-path">Loading...</span></div>
    <ul class="file-list" id="b1-list"></ul>
    <h4>Selected Input Files:</h4>
    <div class="selected-box" id="b1-selected-box">No files selected.</div>
    
    <!-- ASCII Spreadsheet Loader Toggle -->
    <div class="form-group" style="background: #fef3c7; border-color: #fde68a; margin-top: 15px;">
      <label style="cursor: pointer; font-weight: bold;">
        <input type="checkbox" id="toggle-ascii-checkbox" onchange="toggleAsciiLoader(this.checked)" style="width:18px; height:18px; vertical-align:-3px;">
        Enable ASCII Spreadsheet Auto-Loader
      </label>
      <div style="font-size: 13px; color: #4b5563; margin-top: 4px;">Select an ASCII spreadsheet file (*.csv, *.txt, *.tsv) to automatically extract Reflectance, Ancillary, and TOPO Subgroups.</div>
    </div>

    <!-- ASCII Loader Panel -->
    <div id="ascii-loader-container" class="disabled-zone" style="margin-top: 10px; padding: 12px; border: 1px solid #cbd5e1; border-radius: 6px; background: #f8fafc;">
      <h4>Select ASCII Spreadsheet File (*.csv, *.txt, *.tsv, *.dat)</h4>
      <div class="browser-header"><span id="b-ascii-current-path">Loading...</span></div>
      <ul class="file-list" id="b-ascii-list"></ul>
      
      <div id="ascii-preview-section" style="display: none; margin-top: 15px;">
        <h4>ASCII File Preview (First 3 Lines):</h4>
        <pre id="ascii-preview-box" style="background: #1e293b; color: #f1f5f9; padding: 10px; border-radius: 4px; font-size: 12px; overflow-x: auto; font-family: monospace; margin: 0 0 15px 0;"></pre>
        
        <h4>Configure Column Mapping:</h4>
        <div class="grid-3">
          <div>
            <label for="ascii-col-reflectance" style="font-size:12px; font-weight:bold;">Reflectance Column:</label>
            <select id="ascii-col-reflectance" style="width:100%; margin-top:4px;" onchange="applyAsciiData()"></select>
          </div>
          <div>
            <label for="ascii-col-ancillary" style="font-size:12px; font-weight:bold;">Ancillary Column:</label>
            <select id="ascii-col-ancillary" style="width:100%; margin-top:4px;" onchange="applyAsciiData()"></select>
          </div>
          <div>
            <label for="ascii-col-topo" style="font-size:12px; font-weight:bold;">TOPO Subgroup Column:</label>
            <select id="ascii-col-topo" style="width:100%; margin-top:4px;" onchange="applyAsciiData()"></select>
          </div>
        </div>
        <div id="ascii-status-msg" style="margin-top: 10px; font-size: 12px; font-weight: bold; color: #059669;"></div>
      </div>
    </div>
  </div>

  <!-- FILE GROUPING -->
  <div class="card">
    <h3>2. Group Selected Reflectance for TOPO correction subgroups (Optional)</h3>
    <div class="form-group">
      <label for="num-groups"><strong>Number of TOPO subgroups: </strong></label>
      <input type="number" id="num-groups" value="2" min="1" style="width: 80px; margin-right: 10px;" onchange="updateGroupCount(this.value)">
    </div>
    
    <div class="grid-2">
      <div>
        <h4>Define Group Names:</h4>
        <div id="group-names-container" style="display: flex; flex-direction: column; gap: 8px;"></div>
      </div>
      <div>
        <h4>Assign Files to Groups:</h4>
        <div id="file-assignment-container">
          <p style="color: #666; font-size: 14px;">Select files in Step 1 to assign them.</p>
        </div>
      </div>
    </div>
  </div>

  <!-- FILE SELECTOR 3: ANCILLARY IMAGES -->
  <div class="card">
    <h3>3. Ancillary images in ENVI format paired with reflectance images (e.g. *obs_ort, *obs*) Leave it blank for NEON AOP images</h3>
    <div class="browser-header"><span id="b2-current-path">Loading...</span></div>
    <ul class="file-list" id="b2-list"></ul>
    <h4>Selected Ancillary Images:</h4>
    <div class="selected-box" id="b2-selected-box">No files selected.</div>
  </div>

  <!-- SECTION 4: PRECOMPUTED CORRECTION MODEL COEFFICIENTS -->
  <div class="card">
    <h3>4. Precomputed correction model coefficients</h3>
    <div class="form-group" style="background: #eff6ff; border-color: #bfdbfe;">
      <label style="cursor: pointer; font-weight: bold;">
        <input type="checkbox" id="toggle-supp-checkbox" onchange="toggleSupplementary(this.checked)" style="width:18px; height:18px; vertical-align:-3px;">
        Enable Precomputed Correction Models
      </label>
      <div style="font-size: 13px; color: #4b5563; margin-top: 4px;">Check to unlock TOPO and BRDF JSON file selectors below (*.json filter enabled).</div>
    </div>

    <div id="supp-browsers-container" class="disabled-zone">
      <div style="margin-bottom: 20px;">
        <h4>TOPO json files (*.json)</h4>
        <div class="browser-header"><span id="b3-current-path">Loading...</span></div>
        <ul class="file-list" id="b3-list"></ul>
        <h4>Selected TOPO Files:</h4>
        <div class="selected-box" id="b3-selected-box">No files selected.</div>
      </div>

      <div>
        <h4>BRDF json files (*.json)</h4>
        <div class="browser-header"><span id="b4-current-path">Loading...</span></div>
        <ul class="file-list" id="b4-list"></ul>
        <h4>Selected BRDF Files:</h4>
        <div class="selected-box" id="b4-selected-box">No files selected.</div>
      </div>
    </div>
  </div>

  <!-- SECTION 5: USE GLT -->
  <div class="card">
    <h3>5. Use External GLT for location adjustment</h3>
    <div class="form-group" style="background: #f0fdf4; border-color: #bbf7d0;">
      <label style="cursor: pointer; font-weight: bold;">
        <input type="checkbox" id="toggle-glt-checkbox" onchange="toggleGLT(this.checked)" style="width:18px; height:18px; vertical-align:-3px;">
        Enable Geometric Lookup Table (GLT) Processing
      </label>
      <div style="font-size: 13px; color: #4b5563; margin-top: 4px;">Check to unlock the GLT file selector below.</div>
    </div>

    <div id="glt-browser-container" class="disabled-zone">
      <div>
        <h4>Selector GLT (2-band ENVI files paired with reflectance images) </h4>
        <div class="browser-header"><span id="b5-current-path">Loading...</span></div>
        <ul class="file-list" id="b5-list"></ul>
        <h4>Selected GLT Files:</h4>
        <div class="selected-box" id="b5-selected-box">No files selected.</div>
      </div>
    </div>
  </div>

  <!-- SECTION 6: PARAMETERS -->
  <div class="card">
    <h3>6. Processing Parameters</h3>
    <div class="grid-2">
      <!-- Left Column: Radio Groups -->
      <div>
        <div class="form-group">
          <label><strong>Data Type for Input Reflectance Images</strong></label>
          <div class="options-group">
            <label><input type="radio" name="radio1" value="envi" checked> ENVI</label>
            <label><input type="radio" name="radio1" value="neon"> NEON AOP HDF5 (*.h5)</label>
            <label><input type="radio" name="radio1" value="ncav"> AVIRIS NetCDF (*.nc)</label>
            <label><input type="radio" name="radio1" value="emit"> EMIT NetCDF (*.nc)</label>
            <label><input type="radio" name="radio1" value="tanager"> Tanager HDF5 (*.h5)</label>
          </div>
        </div>
        <div class="form-group">
          <label><strong>Data Type for Output Reflectance Images</strong></label>
          <div class="options-group">
            <label><input type="radio" name="radio3" value="envi" checked> ENVI</label>
            <label><input type="radio" name="radio3" value="netcdf"> NetCDF</label>
          </div>
        </div>
        <div class="form-group">
          <label><strong>Good band range in reflectance images</strong></label>
          <div class="options-group">
            <label><input type="radio" name="radio_left_extra" value="badband1" checked> Bands outside ([[300,400],[1337,1430],[1800,1960],[2450,2600]])</label>
            <label><input type="radio" name="radio_left_extra" value="true"> Visible bands (400-750)</label>
            <label><input type="radio" name="radio_left_extra" value="ir"> Bands outside ([[300,750],[1337,1430],[1800,1960],[2450,2600]])</label>
          </div>
        </div>
        <div class="form-group">
          <label><strong>Bands in output reflectance images</strong></label>
          <div class="options-group">
            <label><input type="radio" name="radio2" value="full" checked> Full range</label>
            <label><input type="radio" name="radio2" value="truecolor"> True color ([660,550,440]nm)</label>
            <label><input type="radio" name="radio2" value="subset1"> Selected bands (around [440,550,660,850,976,1650,2217]nm)</label>
          </div>
        </div>        
      </div>
      
      <!-- Right Column: Checkbox Groups -->
      <div>
        <div class="form-group">
          <label><strong>Export Options</strong></label>
          <div class="options-group">
            <label><input type="checkbox" name="check1" value="image" checked> Export Images</label>
            <label><input type="checkbox" name="check1" value="coeffs"> Correction Model Coefficients</label>
            <label><input type="checkbox" name="check1" value="masks"> Save correction related masks</label>
          </div>
        </div>
        <div class="form-group">
          <label><strong>Corrections</strong></label>
          <div class="options-group">
            <label><input type="checkbox" name="check2" value="topo" onchange="updateSuffixDefault()"> Topographic Correction (topo)</label>
            <label><input type="checkbox" name="check2" value="brdf" onchange="updateSuffixDefault()"> BRDF correction (brdf)</label>
            <label><input type="checkbox" name="check2" value="glint" onchange="updateSuffixDefault()"> Sunglint correction (glint)</label>
          </div>
        </div>

      </div>
    </div>

    <!-- THREE DROPDOWNS & ASSOCIATED RADIOS -->
    <div class="grid-3" style="margin-top: 10px;">
      <!-- Dropdown 1 -->
      <div class="form-group">
        <label for="topo_combo"><strong>1. Select TOPO method</strong></label>
        <select id="topo_combo" style="width: 100%; margin-top: 6px;">
          <option value="cosine">Cosine Correction</option>
          <option value="c">C-Correction)</option>
          <option value="scs">Sun-Canopy-Sensor (SCS)</option>
          <option value="scs+c" selected>Sun-Canopy-Sensor with C (SCS+C)</option>
          <option value="mod_minneart">Modified Minnaert</option>
          <option value="precomputed">Precomputed</option>
        </select>
      </div>

      <!-- Dropdown 2 (Modified to contain two radio groups below) -->
      <div class="form-group">
        <label for="brdf_combo"><strong>2. Select BRDF method</strong></label>
        <select id="brdf_combo" style="width: 100%; margin-top: 6px;">
          <option value="universal">Universal BRDF</option>
          <option value="flex" selected>FlexBRDF</option>
          <option value="precomputed">Precomputed</option>
        </select>
        
        <!-- Dropdown 2: Radio Group 1 -->
        <div style="margin-top: 10px; border-top: 1px dashed #d1d5db; padding-top: 8px;">
          <label style="font-size: 12px; color: #374151;"><strong>BRDF grouping mode</strong></label>
          <div class="options-group" style="margin-top: 4px;">
            <label style="font-size: 12px;"><input type="radio" name="radio_drop2_compression" value="group" checked> Same BRDF model for all images</label>
            <label style="font-size: 12px;"><input type="radio" name="radio_drop2_compression" value="single"> Individual model for each image</label>
          </div>
        </div>

        <!-- Dropdown 2: Radio Group 2 (Added) -->
        <div style="margin-top: 10px; border-top: 1px dashed #d1d5db; padding-top: 8px;">
          <label style="font-size: 12px; color: #374151;"><strong>Target solar zenith angle mode</strong></label>
          <div class="options-group" style="margin-top: 4px;">
            <label style="font-size: 12px;"><input type="radio" name="radio_drop2_interleave" value="scene" checked> Normalized to the averaged angle of all images (scene)</label>
            <label style="font-size: 12px;"><input type="radio" name="radio_drop2_interleave" value="line"> Normalized to the averaged angle of each image (line)</label>
          </div>
        </div>
      </div>

      <!-- Dropdown 3 -->
      <div class="form-group">
        <label for="glint_combo"><strong>3. Select Sunglint method</strong></label>
        <select id="glint_combo" style="width: 100%; margin-top: 6px;">
          <option value="hochberg">Hochberg</option>
          <option value="gao" selected>Gao</option>
          <option value="hedley">Hedley</option>
        </select>
        <div style="margin-top: 10px; border-top: 1px dashed #d1d5db; padding-top: 8px;">
          <label style="font-size: 12px; color: #374151;"><strong>Reference wavelength for sunglint correction</strong></label>
          <div class="options-group" style="margin-top: 4px;">
            <label style="font-size: 12px;"><input type="radio" name="radio_drop3" value="860"> 860nm (NIR)</label>
            <label style="font-size: 12px;"><input type="radio" name="radio_drop3" value="1650" checked> 1650nm (SWIR)</label>
            <label style="font-size: 12px;"><input type="radio" name="radio_drop3" value="2190"> 2190nm (SWIR)</label>
          </div>
        </div>
      </div>
    </div>
  </div>

  <!-- EXPORT FOLDER SELECTOR -->
  <div class="card">
    <h3>7. Export Folder Setup (Stored inside JSON Parameters)</h3>
    <p style="font-size: 14px; margin-top:0;">Select the directory where downstream processing scripts should export generated results.</p>
    
    <div class="browser-header">
      <span id="b-export-current-path">Loading...</span>
      <button class="btn-secondary" onclick="setExportDir()" style="flex: initial; padding: 6px 12px; font-size: 13px;">Set as Export Directory</button>
    </div>
    <ul class="file-list" id="b-export-list"></ul>
    <div style="margin-top: 10px;">
      <strong>Selected Export Target Directory:</strong><br>
      <span class="path-highlight" id="selected-export-dir">Not set</span>
    </div>

    <!-- SUFFIX TEXT EDITOR BOX -->
    <div style="margin-top: 15px;">
      <label for="export-suffix"><strong>Output Filename Suffix:</strong></label>
      <div style="font-size: 12px; color: #666; margin-bottom: 4px;">Automatically generated based on chosen Pre-processing Steps in Section 6. Editable below:</div>
      <input type="text" value="raw" id="export-suffix" style="width: 100%;">
    </div>
  </div>

  <!-- JSON FILE SAVE LOCATION SELECTOR -->
  <div class="card">
    <h3>8. Save JSON Configuration File</h3>
    <p style="font-size: 14px; margin-top:0;">Select the server folder and specify the file name for writing the JSON configuration file itself.</p>
    
    <div class="browser-header">
      <span id="b-json-save-current-path">Loading...</span>
      <button class="btn-secondary" onclick="setJsonSaveDir()" style="flex: initial; padding: 6px 12px; font-size: 13px;">Set Save Directory</button>
    </div>
    <ul class="file-list" id="b-json-save-list"></ul>
    
    <div style="margin: 15px 0;" class="grid-2">
      <div>
        <strong>JSON File Save Location:</strong><br>
        <span class="path-highlight" id="selected-json-save-dir">Not set (Defaults to current folder)</span>
      </div>
      <div>
        <label for="json-filename"><strong>Configuration File Name:</strong></label><br>
        <input type="text" id="json-filename" value="server_job_config.json" style="width: 100%; margin-top: 4px;">
      </div>
    </div>

    <!-- ACTION BUTTONS -->
    <div class="btn-group">
      <button onclick="generateJSON()">1. Generate & Process JSON</button>
      <button class="btn-save" onclick="saveJSONToServer()">2. Save JSON File to Server</button>
    </div>
    
    <pre id="json-preview"></pre>
  </div>

  <script>
    const state = {
      b1: { selectedFiles: new Set() },
      b2: { selectedFiles: new Set() },
      b3: { selectedFiles: new Set() },
      b4: { selectedFiles: new Set() },
      b5: { selectedFiles: new Set() }
    };
    
    let currentBrowseDirs = { export: '', jsonSave: '' };
    let selectedExportDir = '';
    let selectedJsonSaveDir = '';
    let isSupplementaryEnabled = false;
    let isGLTEnabled = false;
    let isAsciiLoaderEnabled = false;
    let parsedAsciiData = null;

    function escapeHtml(str) {
      if (!str) return '';
      return String(str)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;")
        .replace(/'/g, "&#039;");
    }
    
    function updateSuffixDefault() {
      const selected = Array.from(document.querySelectorAll('input[name="check2"]:checked')).map(cb => cb.value);
      const allEmpty = selected.every(item => item === "");
      if (allEmpty){
          document.getElementById('export-suffix').value = "raw";
      } else {
          document.getElementById('export-suffix').value = selected.join('_');
      }
    }

    function toggleSupplementary(enabled) {
      isSupplementaryEnabled = enabled;
      const container = document.getElementById('supp-browsers-container');
      container.classList.toggle('disabled-zone', !enabled);
    }

    function toggleGLT(enabled) {
      isGLTEnabled = enabled;
      const container = document.getElementById('glt-browser-container');
      container.classList.toggle('disabled-zone', !enabled);
    }

    function toggleAsciiLoader(enabled) {
      isAsciiLoaderEnabled = enabled;
      const container = document.getElementById('ascii-loader-container');
      container.classList.toggle('disabled-zone', !enabled);
    }

    // Grouping UI Logic
    let numGroups = 2;
    let groupNames = ['Group 1', 'Group 2'];
    let fileAssignments = {}; 

    function updateGroupCount(val) {
      const n = parseInt(val);
      if (isNaN(n) || n < 1) return;
      numGroups = n;
      while (groupNames.length < numGroups) groupNames.push(`Group ${groupNames.length + 1}`);
      if (groupNames.length > numGroups) groupNames = groupNames.slice(0, numGroups);
      renderGroupNamesUI();
      renderFileAssignmentsUI();
    }

    function updateGroupName(index, newName) {
      groupNames[index] = newName;
      renderFileAssignmentsUI();
    }

    function renderGroupNamesUI() {
      const container = document.getElementById('group-names-container');
      container.innerHTML = '';
      groupNames.forEach((name, i) => {
        container.innerHTML += `
          <div style="display: flex; align-items: center; gap: 8px;">
            <label style="width: 25px; font-weight: bold; color: #666;">${i + 1}.</label>
            <input type="text" value="${escapeHtml(name)}" oninput="updateGroupName(${i}, this.value)" style="flex: 1;">
          </div>`;
      });
    }

    function updateFileAssignment(selectEl) {
      const filepath = selectEl.dataset.filepath;
      fileAssignments[filepath] = parseInt(selectEl.value);
    }
    
    function renderFileAssignmentsUI() {
      const container = document.getElementById('file-assignment-container');
      if (state.b1.selectedFiles.size === 0) {
        container.innerHTML = '<p style="color: #666; font-size: 14px; margin:0;">No files selected in Step 1.</p>';
        return;
      }
      let html = '<table class="grouping-table"><thead><tr><th>File Name</th><th>Assigned Group</th></tr></thead><tbody>';
      state.b1.selectedFiles.forEach(filepath => {
        const filename = filepath.replace(/^.*[/\\]/, '');
        if (fileAssignments[filepath] === undefined || fileAssignments[filepath] >= numGroups) fileAssignments[filepath] = 0;
        let optionsHtml = '';
        groupNames.forEach((gname, idx) => {
          optionsHtml += `<option value="${idx}" ${fileAssignments[filepath] === idx ? 'selected' : ''}>${escapeHtml(gname)}</option>`;
        });
        html += `<tr>
            <td title="${escapeHtml(filepath)}" style="max-width: 150px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">📄 ${escapeHtml(filename)}</td>
            <td><select data-filepath="${escapeHtml(filepath)}" onchange="updateFileAssignment(this)" style="width:100%;">${optionsHtml}</select></td>
          </tr>`;
      });
      container.innerHTML = html + '</tbody></table>';
    }

    // Generic Directory Browser (Supports File Extension Filtering)
    async function loadBrowser(path, browserId, selectionType, ext = '') {
      const res = await fetch(`/api/browse?path=${encodeURIComponent(path)}&ext=${encodeURIComponent(ext)}`);
      const data = await res.json();
      
      const currentPathId = `${browserId}-current-path`;
      const listId = `${browserId}-list`;
      
      document.getElementById(currentPathId).textContent = data.current_path;
      
      if (browserId === 'b-export') currentBrowseDirs.export = data.current_path;
      if (browserId === 'b-json-save') currentBrowseDirs.jsonSave = data.current_path;
      
      const list = document.getElementById(listId);
      list.innerHTML = '';

      if (data.parent !== data.current_path) {
        list.innerHTML += `<li><a class="dir-link" data-path="${escapeHtml(data.parent)}" onclick="handleDirClick(this, '${browserId}', '${selectionType}', '${ext}')">📁 .. (Go up)</a></li>`;
      }

      data.dirs.forEach(dir => {
        const dirPath = data.current_path === '/' ? `/${dir}` : `${data.current_path}/${dir}`;
        list.innerHTML += `<li><a class="dir-link" data-path="${escapeHtml(dirPath)}" onclick="handleDirClick(this, '${browserId}', '${selectionType}', '${ext}')">📁 ${escapeHtml(dir)}</a></li>`;
      });

      if (selectionType === 'file') {
        data.files.forEach(file => {
          const filePath = data.current_path === '/' ? `/${file}` : `${data.current_path}/${file}`;
          const isChecked = state[browserId] && state[browserId].selectedFiles.has(filePath) ? 'checked' : '';
          list.innerHTML += `
            <li>
              <label style="cursor:pointer; display:flex; align-items:center; width:100%; user-select:none;">
                <input type="checkbox" style="margin-right:8px;" ${isChecked} data-filepath="${escapeHtml(filePath)}" data-browserid="${browserId}" onchange="handleFileToggle(this)">
                📄 ${escapeHtml(file)}
              </label>
            </li>`;
        });
      } else if (selectionType === 'ascii-file') {
        data.files.forEach(file => {
          const filePath = data.current_path === '/' ? `/${file}` : `${data.current_path}/${file}`;
          list.innerHTML += `
            <li>
              <a class="dir-link" style="color:#059669;" data-path="${escapeHtml(filePath)}" onclick="selectAndParseAscii(this.dataset.path)">📄 ${escapeHtml(file)} (Click to Load)</a>
            </li>`;
        });
      }
    }

    function handleDirClick(el, browserId, selectionType, ext) {
      loadBrowser(el.dataset.path, browserId, selectionType, ext);
    }

    function handleFileToggle(inputEl) {
      toggleFile(inputEl.dataset.filepath, inputEl.checked, inputEl.dataset.browserid);
    }

    async function selectAndParseAscii(filepath) {
      const res = await fetch(`/api/read-ascii?path=${encodeURIComponent(filepath)}`);
      const result = await res.json();
      
      if (!result.success) {
        alert("Error reading file: " + result.message);
        return;
      }

      parsedAsciiData = result;
      
      const previewBox = document.getElementById('ascii-preview-box');
      previewBox.textContent = result.preview.join('\n');
      
      const refSelect = document.getElementById('ascii-col-reflectance');
      const ancSelect = document.getElementById('ascii-col-ancillary');
      const topoSelect = document.getElementById('ascii-col-topo');

      refSelect.innerHTML = '';
      ancSelect.innerHTML = '';
      topoSelect.innerHTML = '';

      result.headers.forEach((hdr, idx) => {
        refSelect.innerHTML += `<option value="${idx}">${escapeHtml(hdr)} (Col ${idx+1})</option>`;
        ancSelect.innerHTML += `<option value="${idx}">${escapeHtml(hdr)} (Col ${idx+1})</option>`;
        topoSelect.innerHTML += `<option value="${idx}">${escapeHtml(hdr)} (Col ${idx+1})</option>`;
      });

      result.headers.forEach((hdr, idx) => {
        const lowerHdr = hdr.toLowerCase();
        if (lowerHdr.includes('ref') || lowerHdr.includes('primary') || lowerHdr.includes('image')) refSelect.value = idx;
        if (lowerHdr.includes('anc') || lowerHdr.includes('obs') || lowerHdr.includes('obs_ort')) ancSelect.value = idx;
        if (lowerHdr.includes('topo') || lowerHdr.includes('group') || lowerHdr.includes('sub')) topoSelect.value = idx;
      });

      document.getElementById('ascii-preview-section').style.display = 'block';
      applyAsciiData();
    }

    function applyAsciiData() {
      if (!parsedAsciiData || !parsedAsciiData.rows) return;

      const refColIdx = parseInt(document.getElementById('ascii-col-reflectance').value);
      const ancColIdx = parseInt(document.getElementById('ascii-col-ancillary').value);
      const topoColIdx = parseInt(document.getElementById('ascii-col-topo').value);

      if (isNaN(refColIdx) || isNaN(ancColIdx) || isNaN(topoColIdx)) return;

      state.b1.selectedFiles.clear();
      parsedAsciiData.rows.forEach(row => {
        if (row[refColIdx]) state.b1.selectedFiles.add(row[refColIdx]);
      });
      document.getElementById('b1-selected-box').innerHTML = state.b1.selectedFiles.size === 0 ? 'No files selected.' : Array.from(state.b1.selectedFiles).map(escapeHtml).join('<br>');

      state.b2.selectedFiles.clear();
      parsedAsciiData.rows.forEach(row => {
        if (row[ancColIdx]) state.b2.selectedFiles.add(row[ancColIdx]);
      });
      document.getElementById('b2-selected-box').innerHTML = state.b2.selectedFiles.size === 0 ? 'No files selected.' : Array.from(state.b2.selectedFiles).map(escapeHtml).join('<br>');

      const uniqueSubgroups = [];
      const fileToSubgroupMap = {};

      parsedAsciiData.rows.forEach(row => {
        const refFile = row[refColIdx];
        const subgroup = row[topoColIdx] || 'Unassigned';
        if (refFile) {
          if (!uniqueSubgroups.includes(subgroup)) {
            uniqueSubgroups.push(subgroup);
          }
          fileToSubgroupMap[refFile] = subgroup;
        }
      });

      if (uniqueSubgroups.length > 0) {
        groupNames = uniqueSubgroups;
        numGroups = uniqueSubgroups.length;
        document.getElementById('num-groups').value = numGroups;

        fileAssignments = {};
        Object.entries(fileToSubgroupMap).forEach(([refFile, subgroup]) => {
          fileAssignments[refFile] = groupNames.indexOf(subgroup);
        });

        renderGroupNamesUI();
        renderFileAssignmentsUI();
      }

      document.getElementById('ascii-status-msg').textContent = `✓ Loaded ${state.b1.selectedFiles.size} Reflectance files, ${state.b2.selectedFiles.size} Ancillary files, and ${uniqueSubgroups.length} TOPO Subgroups.`;
    }

    function toggleFile(path, isChecked, browserId) {
      if (isChecked) state[browserId].selectedFiles.add(path);
      else state[browserId].selectedFiles.delete(path);
      
      const box = document.getElementById(`${browserId}-selected-box`);
      box.innerHTML = state[browserId].selectedFiles.size === 0 ? 'No files selected.' : Array.from(state[browserId].selectedFiles).map(escapeHtml).join('<br>');
      if (browserId === 'b1') renderFileAssignmentsUI();
    }

    function setExportDir() {
      selectedExportDir = currentBrowseDirs.export;
      document.getElementById('selected-export-dir').textContent = selectedExportDir;
    }

    function setJsonSaveDir() {
      selectedJsonSaveDir = currentBrowseDirs.jsonSave;
      document.getElementById('selected-json-save-dir').textContent = selectedJsonSaveDir;
    }

    function collectRawPayload() {
      const groupedFiles = {};
      groupNames.forEach(name => { groupedFiles[name] = []; });
      state.b1.selectedFiles.forEach(filepath => {
        const groupIdx = fileAssignments[filepath];
        if (groupIdx !== undefined && groupNames[groupIdx]) {
          groupedFiles[groupNames[groupIdx]].push(filepath);
        }
      });

      const selectedTOPOMethod = document.getElementById('topo_combo').value;
      const selectedBRDFMethod = document.getElementById('brdf_combo').value;
      const selectedGlintMethod = document.getElementById('glint_combo').value;
      const exportSuffix = document.getElementById('export-suffix').value;

      return {
        timestamp: new Date().toISOString(),
        export_directory: selectedExportDir,        
        inputs: {
          images: Array.from(state.b1.selectedFiles),
          primary_files_grouped: groupedFiles,
          anc_files: Array.from(state.b2.selectedFiles),
          enable_precomputed_correction: isSupplementaryEnabled,
          enable_ascii_loader: isAsciiLoaderEnabled,
          topo_json_files: isSupplementaryEnabled ? Array.from(state.b3.selectedFiles) : [],
          brdf_json_files: isSupplementaryEnabled ? Array.from(state.b4.selectedFiles) : [],
          use_glt: isGLTEnabled,
          glt_files: isGLTEnabled ? Array.from(state.b5.selectedFiles) : []
        },
        parameters: {
          data_type: document.querySelector('input[name="radio1"]:checked').value,
          data_type_export: document.querySelector('input[name="radio3"]:checked').value,
          subset_waves: document.querySelector('input[name="radio2"]:checked').value,
          bad_bands: document.querySelector('input[name="radio_left_extra"]:checked').value,
          tasks: Array.from(document.querySelectorAll('input[name="check1"]:checked')).map(cb => cb.value),
          corrections: Array.from(document.querySelectorAll('input[name="check2"]:checked')).map(cb => cb.value),
          suffix: exportSuffix,
          topo_method: selectedTOPOMethod,
          brdf_method: selectedBRDFMethod,
          brdf_group_mode: document.querySelector('input[name="radio_drop2_compression"]:checked').value,
          brdf_solar_zn: document.querySelector('input[name="radio_drop2_interleave"]:checked').value,
          glint_method: selectedGlintMethod,
          glint_ref_wave: document.querySelector('input[name="radio_drop3"]:checked').value
        }
      };
    }

    // 1. Generate & Process Preview
    async function generateJSON() {
      const rawPayload = collectRawPayload();
      
      const res = await fetch('/api/process-preview', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(rawPayload)
      });
      
      const processedConfig = await res.json();
      const previewEl = document.getElementById('json-preview');
      previewEl.style.display = 'block';
      previewEl.style.color = '#38bdf8';
      previewEl.textContent = `// Processed & Reordered JSON Config (Internal Python Transformation):\n\n` + JSON.stringify(processedConfig, null, 2);
      window.scrollTo({ top: document.body.scrollHeight, behavior: 'smooth' });
    }

    // 2. Save JSON to File System
    async function saveJSONToServer() {
      const rawPayload = collectRawPayload();
      const filename = document.getElementById('json-filename').value.trim() || 'server_job_config.json';

      const requestData = {
        raw_payload: rawPayload,
        save_directory: selectedJsonSaveDir,
        filename: filename
      };

      try {
        const res = await fetch('/api/save', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(requestData)
        });
        
        const result = await res.json();
        const previewEl = document.getElementById('json-preview');
        previewEl.style.display = 'block';
        
        if (result.status === 'success') {
          previewEl.style.color = '#a7f3d0';
          previewEl.textContent = `// SUCCESS: File written to disk!\n// JSON Save Path: ${result.saved_file_path}\n\n` + JSON.stringify(result.saved_config, null, 2);
        } else {
          previewEl.style.color = '#fca5a5';
          previewEl.textContent = `// ERROR: ${result.message}`;
        }
        window.scrollTo({ top: document.body.scrollHeight, behavior: 'smooth' });
      } catch (err) {
        alert("Failed to communicate with server.");
      }
    }

    // Initialize UI
    renderGroupNamesUI();
    renderFileAssignmentsUI();
    updateSuffixDefault();
    loadBrowser('', 'b1', 'file');
    loadBrowser('', 'b-ascii', 'ascii-file');
    loadBrowser('', 'b2', 'file');
    loadBrowser('', 'b3', 'file', '.json');
    loadBrowser('', 'b4', 'file', '.json');
    loadBrowser('', 'b5', 'file');
    loadBrowser('', 'b-export', 'dir');
    loadBrowser('', 'b-json-save', 'dir');
  </script>
</body>
</html>
"""

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/api/browse')
def browse():
    req_path = request.args.get('path')
    ext_filter = request.args.get('ext', '')

    if not req_path:
        req_path = os.path.expanduser('~')

    target_path = os.path.abspath(req_path)
    parent_path = os.path.dirname(target_path)
    dirs, files = [], []

    try:
        for entry in os.scandir(target_path):
            if entry.is_dir():
                dirs.append(entry.name)
            elif entry.is_file():
                if ext_filter:
                    if entry.name.lower().endswith(ext_filter.lower()):
                        files.append(entry.name)
                else:
                    files.append(entry.name)
    except PermissionError:
        pass

    dirs.sort()
    files.sort()
    return jsonify({'current_path': target_path, 'parent': parent_path, 'dirs': dirs, 'files': files})

@app.route('/api/read-ascii')
def read_ascii():
    filepath = request.args.get('path')
    if not filepath or not os.path.isfile(filepath):
        return jsonify({'success': False, 'message': 'File not found'}), 400

    try:
        with open(filepath, 'r', encoding='utf-8', errors='strict') as f:
            lines = []
            for _ in range(50):
                line = f.readline()
                if not line:
                    break
                lines.append(line)

        raw_preview = [l.strip('\r\n') for l in lines[:3] if l.strip()]
        if not raw_preview:
            return jsonify({'success': False, 'message': 'File is empty'}), 400

        header_line = raw_preview[0]
        if ',' in header_line:
            delimiter = ','
        elif '\t' in header_line:
            delimiter = '\t'
        elif ';' in header_line:
            delimiter = ';'
        else:
            delimiter = None

        with open(filepath, 'r', encoding='utf-8') as f:
            if delimiter:
                reader = csv.reader(f, delimiter=delimiter)
            else:
                reader = csv.reader(f, delimiter=' ', skipinitialspace=True)

            all_rows = [row for row in reader if row and not row[0].startswith('#')]

        if not all_rows:
            return jsonify({'success': False, 'message': 'No valid data rows found in ASCII file'}), 400

        headers = [h.strip() for h in all_rows[0]]
        data_rows = [[cell.strip() for cell in row] for row in all_rows[1:]]

        return jsonify({
            'success': True,
            'preview': raw_preview,
            'headers': headers,
            'rows': data_rows
        })

    except UnicodeDecodeError:
        return jsonify({'success': False, 'message': 'Not a valid ASCII text file'}), 400
    except Exception as e:
        return jsonify({'success': False, 'message': str(e)}), 500


@app.route('/api/process-preview', methods=['POST'])
def process_preview():
    raw_payload = request.json
    processed_config = process_and_reorder_config(raw_payload)
    return jsonify(processed_config)

@app.route('/api/save', methods=['POST'])
def save():
    data = request.json

    raw_payload = data.get('raw_payload', {})
    save_dir = data.get('save_directory', '')
    filename = data.get('filename', 'server_job_config.json')

    if not filename.endswith('.json'):
        filename += '.json'

    final_config = process_and_reorder_config(raw_payload)

    if len(final_config["job_metadata"]["warning_message"])>0:
        return jsonify({"status": "Warning", "message": "\n".join(final_config["job_metadata"]["warning_message"])}), 500

    if not save_dir or not os.path.isdir(save_dir):
        save_dir = os.path.abspath(os.path.dirname(__file__))

    full_save_path = os.path.join(save_dir, filename)

    try:
        with open(full_save_path, 'w') as f:
            json.dump(final_config, f, indent=4)
        return jsonify({
            "status": "success",
            "saved_file_path": full_save_path,
            "saved_config": final_config
        })
    except Exception as e:
        return jsonify({"status": "error", "message": str(e)}), 500

if __name__ == '__main__':
    if len(sys.argv)>1:
        port_user = int(sys.argv[1])
    elif len(sys.argv)==1:
        port_user = 5005
    app.run(host='127.0.0.1', port=port_user, debug=True)

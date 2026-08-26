import os
import sys
import json
import csv
from datetime import datetime, timezone
import numpy as np
from flask import Flask, request, jsonify

# ==============================================================================
# CORE CONFIGURATION LOGIC (Ported directly from Tkinter source)
# ==============================================================================

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

    if len(images) < 1:
        feedback_message += ['// No reflectance images is selected.']

    aviris_anc_names = ['path_length','sensor_az','sensor_zn',
                        'solar_az', 'solar_zn','phase','slope',
                        'aspect', 'cosine_i','utc_time']

    if config_dict["file_type"] in ["envi","emit"]:
        config_dict["anc_files"] = {}

        if ascii_loader_bool:
            anc_files = setting_dict["input_datasets"]["anc_files"]
        else:    
            anc_files = sorted(setting_dict["input_datasets"]["anc_files"])

        if 'topo' in corr_list or 'brdf' in corr_list:
            if not len(images) == len(anc_files):
                feedback_message += ['// Reflectance images and ancillary images do not match in number.']
            elif len(anc_files) == 0:
                feedback_message += ['// Ancillary images required for TOPO or BRDF correction, but none is provided.']
            else:
                for i, image in enumerate(images):
                    config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                                [[anc_files[i], a] for a in range(len(aviris_anc_names))]))

    elif config_dict["file_type"] == "neon":
        pass
    elif config_dict["file_type"] == "ncav":
        ncav_anc_names = ['path_length','to_sensor_azimuth','to_sensor_zenith',
                        'to_sun_azimuth', 'to_sun_zenith','solar_phase','slope',
                        'aspect', 'cosine_i','utc_time']

        config_dict["anc_files"] = {}

        if ascii_loader_bool:
            anc_files = setting_dict["input_datasets"]["anc_files"]
        else:    
            anc_files = sorted(setting_dict["input_datasets"]["anc_files"])

        if 'topo' in corr_list or 'brdf' in corr_list:
            if not len(images) == len(anc_files):
                feedback_message += ['// Reflectance images and ancillary images do not match in number.']
            elif len(anc_files) == 0:
                feedback_message += ['// Ancillary images required for TOPO or BRDF correction, but none is provided.']
            else:
                for i, image in enumerate(images):
                    config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                                [[anc_files[i], ncav_anc_names[a]] for a in range(len(aviris_anc_names))]))

    config_dict['output_dir'] = setting_dict["export_settings"]["target_export_directory"] + "/"
    config_dict['export_type'] = setting_dict["processing_pipeline"]["export_type"]
    config_dict["use_glt"] = setting_dict["input_datasets"]["use_glt"]

    if config_dict["use_glt"]:
        if ascii_loader_bool:
            glt_files = setting_dict["input_datasets"]["glt_files"]
        else:
            glt_files = sorted(setting_dict["input_datasets"]["glt_files"])
        config_dict["glt_files"] = {}
        if not len(images) == len(glt_files):
            feedback_message += ['// Reflectance images and GLT files do not match.']
        elif len(glt_files) == 0:
            feedback_message += ['// GLT files required, but none is provided.']
        else:
            for i, image in enumerate(images):
                config_dict['glt_files'][image] = dict(zip(["glt_x","glt_y"],
                                                [[glt_files[i], a] for a in range(2)]))

    config_dict["corrections"] = corr_list
    config_dict["chunk_size"] = setting_dict["processing_pipeline"]["chunk_size"]

    config_dict["topo"] = {}
    if ascii_loader_bool:
        topo_json_files = setting_dict["input_datasets"]["precomputed"]["topo_json_files"]
    else:
        topo_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["topo_json_files"])

    if 'topo' in corr_list:
        if not len(images) == len(topo_json_files):
            feedback_message += ['// Reflectance images and TOPO coefficient files do not match.']
        elif len(topo_json_files) == 0:
            feedback_message += ['// Precomputed TOPO coefficient files required, but none is provided.']
        else:
            config_dict["topo"] = dict(zip(images, topo_json_files))

    config_dict["brdf"] = {}
    if ascii_loader_bool:
        brdf_json_files = setting_dict["input_datasets"]["precomputed"]["brdf_json_files"]
    else:
        brdf_json_files = sorted(setting_dict["input_datasets"]["precomputed"]["brdf_json_files"])

    if 'brdf' in corr_list:
        if not len(images) == len(brdf_json_files):
            feedback_message += ['// Reflectance images and BRDF coefficient files do not match.']
        elif len(brdf_json_files) == 0:
            feedback_message += ['// Precomputed BRDF coefficient files required, but none is provided.']
        else:
            config_dict["brdf"] = dict(zip(images, brdf_json_files))

    config_dict["resampling"] = {}
    if not setting_dict["processing_pipeline"]["resample_method"] == 'None':        
        config_dict["resampling"]['type'] = setting_dict["processing_pipeline"]["resample_method"]

    config_dict["masks"] = [["ndi", {'band_1': 850,'band_2': 660, 'min': 0.1,'max': 1.0}]]

    if config_dict["file_type"] == "neon":
        config_dict["masks"] += [['neon_edge', {'radius': 30}]]                        

    config_dict["trait_models"] = setting_dict["input_datasets"]["trait_files"]

    return config_dict, feedback_message


def process_and_reorder_config(raw_data):
    inputs = raw_data.get('inputs', {})
    params = raw_data.get('parameters', {})

    bad_band_range = {
        "badband1": [[300, 400], [1337, 1430], [1800, 1960], [2450, 2600]],
        "true": [[300, 400], [750, 2600]],
        "ir": [[300, 750], [1337, 1430], [1800, 1960], [2450, 2600]]
    }

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
            "chunk_size": params.get('chunk_size', []),
            "bad_bands": bad_band_range.get(params.get('bad_bands'), []),
            "resample_method": params.get('resample_method')
        }
    }
    structured_config, feedback_message = generate_config(structured_config_from_gui)
    return structured_config, feedback_message

# ==============================================================================
# FLASK WEB APP & API ROUTING
# ==============================================================================

app = Flask(__name__)

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Remote Configurator for Trait Mapping</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <style>
        .custom-scrollbar::-webkit-scrollbar { width: 8px; height: 8px; }
        .custom-scrollbar::-webkit-scrollbar-track { background: #1f2937; }
        .custom-scrollbar::-webkit-scrollbar-thumb { background: #4b5563; border-radius: 4px; }
        .list-container { min-height: 80px; max-height: 150px; overflow-y: auto; background: #1f2937; color: #a7f3d0; padding: 0.5rem; font-family: monospace; }
        .hidden { display: none !important; }
    </style>
</head>
<body class="bg-gray-100 text-gray-800 font-sans p-6 pb-20">
    <div class="max-w-5xl mx-auto space-y-6">
        <h1 class="text-2xl font-bold mb-4 text-center">Remote Server Job Configurator for Trait Mapping by HyTools</h1>

        <!-- 1. Primary Files -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">1. Input Reflectance Images (Required)</h2>
            <div class="mb-3">
                <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded" onclick="openFsModal('primary_files', 'multi')">Browse Reflectance Files...</button>
                <button type="button" class="bg-gray-400 hover:bg-gray-500 text-white px-3 py-1 rounded ml-2" onclick="clearList('primary_files')">Clear Files</button>
            </div>
            <ul id="primary_files" class="list-container rounded custom-scrollbar"></ul>

            <!-- ASCII Loader -->
            <div class="mt-4 p-4 border border-gray-300 bg-gray-50 rounded">
                <label class="flex items-center space-x-2 font-semibold">
                    <input type="checkbox" id="enable_ascii" onchange="toggleAscii(this)">
                    <span>Enable ASCII Spreadsheet Auto-Loader (*.csv, *.txt, *.tsv)</span>
                </label>
                <div id="ascii_section" class="mt-3 hidden space-y-3">
                    <button type="button" class="bg-indigo-600 hover:bg-indigo-700 text-white px-3 py-1 rounded" onclick="openFsModal('ascii_file', 'single')">Browse ASCII File...</button>
                    <input type="hidden" id="ascii_file">
                    <span id="ascii_file_display" class="ml-2 text-sm text-gray-600 font-mono"></span>
                    
                    <div>
                        <h3 class="text-sm font-bold">File Preview (First 3 Lines):</h3>
                        <pre id="ascii_preview" class="bg-gray-800 text-gray-100 p-2 text-xs rounded mt-1 overflow-x-auto"></pre>
                    </div>

                    <div class="grid grid-cols-2 gap-4 text-sm mt-2" id="ascii_mapping_ui">
                        <div><label class="block font-bold">Reflectance Column:</label><select id="col_ref" class="w-full border p-1"></select></div>
                        <div><label class="block font-bold">Ancillary Column:</label><select id="col_anc" class="w-full border p-1"></select></div>
                        <div><label class="block font-bold">GLT Column:</label><select id="col_glt" class="w-full border p-1"></select></div>
                        <div><label class="block font-bold">TOPO Coeff Column:</label><select id="col_topo" class="w-full border p-1"></select></div>
                        <div><label class="block font-bold">BRDF Coeff Column:</label><select id="col_brdf" class="w-full border p-1"></select></div>
                    </div>
                    <button type="button" class="bg-green-600 hover:bg-green-700 text-white px-4 py-1 rounded mt-2" onclick="applyAscii()">Load ASCII Data</button>
                    <p id="ascii_status" class="text-green-600 font-bold text-sm"></p>
                </div>
            </div>
        </div>

        <!-- 2. Ancillary Data -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">2. (Optional) Ancillary images in ENVI format</h2>
            <div class="mb-3">
                <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded" onclick="openFsModal('anc_files', 'multi')">Browse Ancillary Files...</button>
                <button type="button" class="bg-gray-400 hover:bg-gray-500 text-white px-3 py-1 rounded ml-2" onclick="clearList('anc_files')">Clear Files</button>
            </div>
            <ul id="anc_files" class="list-container rounded custom-scrollbar"></ul>
        </div>

        <!-- 3. Precomputed Correction -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">3. Precomputed correction model coefficients</h2>
            <label class="flex items-center space-x-2 font-semibold mb-3">
                <input type="checkbox" id="enable_precomputed" onchange="togglePrecomputed(this)">
                <span>Use precomputed correction model coefficients</span>
            </label>
            <div id="precomputed_section" class="hidden space-y-4">
                <div>
                    <label class="block text-sm font-bold mb-1">TOPO Coefficients (JSON)</label>
                    <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded mb-2 text-sm" onclick="openFsModal('topo_files', 'multi')">Browse TOPO JSON...</button>
                    <button type="button" class="bg-gray-400 text-white px-3 py-1 rounded mb-2 text-sm" onclick="clearList('topo_files')">Clear</button>
                    <ul id="topo_files" class="list-container rounded custom-scrollbar"></ul>
                </div>
                <div>
                    <label class="block text-sm font-bold mb-1">BRDF Coefficients (JSON)</label>
                    <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded mb-2 text-sm" onclick="openFsModal('brdf_files', 'multi')">Browse BRDF JSON...</button>
                    <button type="button" class="bg-gray-400 text-white px-3 py-1 rounded mb-2 text-sm" onclick="clearList('brdf_files')">Clear</button>
                    <ul id="brdf_files" class="list-container rounded custom-scrollbar"></ul>
                </div>
            </div>
        </div>

        <!-- 4. GLT -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">4. Use External GLT for location adjustment</h2>
            <label class="flex items-center space-x-2 font-semibold mb-3">
                <input type="checkbox" id="enable_glt" onchange="toggleGlt(this)">
                <span>Use GLT file</span>
            </label>
            <div id="glt_section" class="hidden">
                <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded mb-2 text-sm" onclick="openFsModal('glt_files', 'multi')">Browse GLT Files...</button>
                <button type="button" class="bg-gray-400 text-white px-3 py-1 rounded mb-2 text-sm" onclick="clearList('glt_files')">Clear</button>
                <ul id="glt_files" class="list-container rounded custom-scrollbar"></ul>
            </div>
        </div>

        <!-- 5. Settings -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">5. Correction Settings</h2>
            <div class="grid grid-cols-2 gap-6">
                <!-- Col 1 -->
                <div>
                    <h3 class="font-bold text-sm mb-1">Data Type for Input Reflectance Images</h3>
                    <select id="data_type" class="w-full border p-2 rounded mb-4">
                        <option value="envi" selected>ENVI</option>
                        <option value="neon">NEON AOP HDF5 (*.h5)</option>
                        <option value="ncav">AVIRIS NetCDF (*.nc)</option>
                        <option value="emit">EMIT NetCDF (*.nc)</option>
                        <option value="tanager">Tanager HDF5 (*.h5)</option>
                    </select>

                    <h3 class="font-bold text-sm mb-1">Data Type for Output Reflectance Images</h3>
                    <select id="data_type_export" class="w-full border p-2 rounded mb-4">
                        <option value="envi" selected>ENVI</option>
                        <option value="netcdf">NetCDF</option>
                    </select>

                    <h3 class="font-bold text-sm mb-1">Good band range in reflectance images</h3>
                    <select id="bad_bands" class="w-full border p-2 rounded mb-4">
                        <option value="badband1" selected>Bands outside ([[300,400],[1337,1430],[1800,1960],[2450,2600]])</option>
                        <option value="true">Visible bands (400-750)</option>
                        <option value="ir">Bands outside ([[300,750],[1337,1430],[1800,1960],[2450,2600]])</option>
                    </select>
                </div>
                <!-- Col 2 -->
                <div>
                    <h3 class="font-bold text-sm mb-1">Correction</h3>
                    <div class="mb-4">
                        <label class="flex items-center space-x-2"><input type="checkbox" id="corr_topo" value="topo"> <span>Topographic Correction (topo)</span></label>
                        <label class="flex items-center space-x-2"><input type="checkbox" id="corr_brdf" value="brdf"> <span>BRDF correction (brdf)</span></label>
                    </div>

                    <div class="flex space-x-4 mb-4">
                        <div class="w-1/2">
                            <h3 class="font-bold text-sm mb-1">Chunk Row Size</h3>
                            <select id="chunk_row" class="w-full border p-2 rounded">
                                <option value="row">Full Row</option>
                                <option value="16">16</option>
                                <option value="64" selected>64</option>
                                <option value="128">128</option>
                                <option value="256">256</option>
                            </select>
                        </div>
                        <div class="w-1/2">
                            <h3 class="font-bold text-sm mb-1">Chunk Col Size</h3>
                            <select id="chunk_col" class="w-full border p-2 rounded">
                                <option value="col">Full Column</option>
                                <option value="16">16</option>
                                <option value="64" selected>64</option>
                                <option value="128">128</option>
                                <option value="256">256</option>
                            </select>
                        </div>
                    </div>

                    <h3 class="font-bold text-sm mb-1">Select Resampling Method</h3>
                    <select id="resample_method" class="w-full border p-2 rounded mb-4">
                        <option value="None">None</option>
                        <option value="nearest" selected>nearest</option>
                        <option value="linear">linear</option>
                        <option value="cubic">cubic</option>
                    </select>
                </div>
            </div>
        </div>

        <!-- 6. Trait Models -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">6. Select Trait Model JSON File</h2>
            <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded mb-2 text-sm" onclick="openFsModal('trait_files', 'multi')">Browse Trait JSON Files...</button>
            <button type="button" class="bg-gray-400 text-white px-3 py-1 rounded mb-2 text-sm" onclick="clearList('trait_files')">Clear</button>
            <ul id="trait_files" class="list-container rounded custom-scrollbar"></ul>
        </div>

        <!-- 7. Export Settings -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">7. Export Folder Setup</h2>
            <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded mb-2" onclick="openFsModal('export_dir', 'dir')">Select Export Directory...</button>
            <input type="text" id="export_dir" readonly class="w-full p-2 border font-mono text-fuchsia-600 bg-gray-50" value="Not set">
        </div>

        <!-- 8. Save JSON -->
        <div class="bg-white p-5 rounded shadow">
            <h2 class="text-lg font-semibold mb-3 border-b pb-2">8. Save JSON Configuration File</h2>
            <div class="grid grid-cols-2 gap-4 mb-4">
                <div>
                    <label class="block font-bold mb-1">JSON Save Folder:</label>
                    <button type="button" class="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded mb-2" onclick="openFsModal('save_dir', 'dir')">Select JSON Save Folder...</button>
                    <input type="text" id="save_dir" readonly class="w-full p-2 border font-mono text-fuchsia-600 bg-gray-50">
                </div>
                <div>
                    <label class="block font-bold mb-1">Configuration File Name:</label>
                    <input type="text" id="save_filename" class="w-full p-2 border mt-8" value="local_job_config.json">
                </div>
            </div>

            <div class="flex space-x-4">
                <button type="button" class="flex-1 bg-gray-800 hover:bg-black text-white font-bold py-3 rounded" onclick="generatePreview()">1. Generate & Preview JSON</button>
                <button type="button" class="flex-1 bg-green-600 hover:bg-green-700 text-white font-bold py-3 rounded" onclick="saveJson()">2. Save JSON to Disk</button>
            </div>
            
            <pre id="json_preview" class="mt-4 bg-gray-900 text-sky-400 p-4 rounded min-h-[200px] overflow-auto custom-scrollbar text-sm font-mono"></pre>
        </div>
    </div>

    <!-- Server File Browser Modal -->
    <div id="fsModal" class="hidden fixed inset-0 bg-gray-900 bg-opacity-75 flex items-center justify-center z-50">
        <div class="bg-white w-3/4 max-w-4xl h-[80vh] flex flex-col rounded-lg shadow-xl">
            <div class="p-4 border-b flex justify-between items-center bg-gray-100 rounded-t-lg">
                <h3 class="text-xl font-bold" id="fsModalTitle">Server File Browser</h3>
                <button class="text-red-500 font-bold text-2xl leading-none" onclick="closeFsModal()">&times;</button>
            </div>
            <div class="p-4 bg-gray-50 flex items-center border-b">
                <button class="bg-gray-300 hover:bg-gray-400 px-3 py-1 rounded mr-3" onclick="loadFsParent()">&#8593; Up</button>
                <input type="text" id="fsCurrentPath" class="flex-1 border p-1 px-2 font-mono" readonly>
            </div>
            <div class="flex-1 overflow-y-auto p-4 custom-scrollbar" id="fsItemList"></div>
            <div class="p-4 border-t bg-gray-100 flex justify-end space-x-3 rounded-b-lg">
                <button class="bg-gray-400 text-white px-4 py-2 rounded" onclick="closeFsModal()">Cancel</button>
                <button class="bg-blue-600 text-white px-4 py-2 rounded font-bold" onclick="confirmFsSelection()">Confirm Selection</button>
            </div>
        </div>
    </div>

    <script>
        // State
        let fsTarget = null;
        let fsMode = 'multi'; // 'multi', 'single', 'dir'
        let currentPath = '';
        let selectedItems = new Set();
        let listData = {
            'primary_files': [], 'anc_files': [], 'topo_files': [], 
            'brdf_files': [], 'glt_files': [], 'trait_files': []
        };

        // --- File System Browser UI ---
        function openFsModal(target, mode) {
            fsTarget = target;
            fsMode = mode;
            selectedItems.clear();
            document.getElementById('fsModal').classList.remove('hidden');
            
            let startPath = '/';
            // Context sensitive starting paths
            if(fsMode === 'dir' && document.getElementById(target).value && document.getElementById(target).value !== 'Not set') {
                startPath = document.getElementById(target).value;
            } else if (target === 'save_dir') {
                startPath = '.'; // backend converts to cwd
            }
            loadFs(startPath);
        }

        function closeFsModal() {
            document.getElementById('fsModal').classList.add('hidden');
        }

        function loadFs(path) {
            fetch(`/api/fs?path=${encodeURIComponent(path)}`)
                .then(r => r.json())
                .then(data => {
                    if(data.error) { alert('Error: ' + data.error); return; }
                    currentPath = data.path;
                    document.getElementById('fsCurrentPath').value = currentPath;
                    renderFsItems(data.dirs, data.files);
                }).catch(e => alert("Failed to load path"));
        }

        function loadFsParent() {
            // Very naive parent calculation, better to rely on API if implemented. 
            // Here backend /api/fs handles the resolution
            let parts = currentPath.replace(/\/$/, '').split('/');
            parts.pop();
            let parent = parts.join('/') || '/';
            loadFs(parent);
        }

        function renderFsItems(dirs, files) {
            const list = document.getElementById('fsItemList');
            list.innerHTML = '';
            
            // Directories
            dirs.forEach(d => {
                const div = document.createElement('div');
                div.className = "p-2 hover:bg-gray-200 cursor-pointer flex items-center";
                div.innerHTML = `<span class="mr-2 text-yellow-500">📁</span> ${d.name}`;
                div.onclick = () => {
                    if (fsMode === 'dir') {
                        if(selectedItems.has(d.path)) {
                            selectedItems.delete(d.path);
                            div.classList.remove('bg-blue-100');
                        } else {
                            selectedItems.clear(); // only one dir selected
                            selectedItems.add(d.path);
                            renderFsItems(dirs, files); // re-render to clear other highlights
                        }
                    } else {
                        loadFs(d.path); // navigate
                    }
                };
                if(fsMode === 'dir' && selectedItems.has(d.path)) div.classList.add('bg-blue-100', 'font-bold');
                list.appendChild(div);
            });

            // Files
            if(fsMode !== 'dir') {
                files.forEach(f => {
                    const div = document.createElement('div');
                    div.className = `p-2 cursor-pointer flex items-center border-t border-gray-100 ${selectedItems.has(f.path) ? 'bg-blue-100 font-bold' : 'hover:bg-gray-100'}`;
                    div.innerHTML = `<span class="mr-2 text-gray-500">📄</span> ${f.name}`;
                    div.onclick = () => {
                        if(fsMode === 'single') {
                            selectedItems.clear();
                            selectedItems.add(f.path);
                            renderFsItems(dirs, files);
                        } else {
                            if(selectedItems.has(f.path)) selectedItems.delete(f.path);
                            else selectedItems.add(f.path);
                            div.classList.toggle('bg-blue-100');
                            div.classList.toggle('font-bold');
                        }
                    };
                    list.appendChild(div);
                });
            }
        }

        function confirmFsSelection() {
            const items = Array.from(selectedItems);
            if (fsMode === 'dir' || fsMode === 'single') {
                if(items.length > 0) {
                    if (fsMode === 'single' && fsTarget === 'ascii_file') {
                        document.getElementById('ascii_file').value = items[0];
                        document.getElementById('ascii_file_display').innerText = items[0];
                        parseAsciiPreview(items[0]);
                    } else {
                        document.getElementById(fsTarget).value = items[0];
                    }
                }
            } else { // Multi file arrays
                items.forEach(item => {
                    if(!listData[fsTarget].includes(item)) listData[fsTarget].push(item);
                });
                updateListUI(fsTarget);
            }
            closeFsModal();
        }

        function updateListUI(targetId) {
            const ul = document.getElementById(targetId);
            ul.innerHTML = '';
            listData[targetId].forEach(item => {
                const li = document.createElement('li');
                li.innerText = item;
                li.className = 'whitespace-nowrap';
                ul.appendChild(li);
            });
        }

        function clearList(targetId) {
            listData[targetId] = [];
            updateListUI(targetId);
        }

        // --- View Toggles ---
        function toggleAscii(cb) { document.getElementById('ascii_section').classList.toggle('hidden', !cb.checked); }
        function togglePrecomputed(cb) { document.getElementById('precomputed_section').classList.toggle('hidden', !cb.checked); }
        function toggleGlt(cb) { document.getElementById('glt_section').classList.toggle('hidden', !cb.checked); }

        // --- ASCII Processing ---
        function parseAsciiPreview(filepath) {
            fetch('/api/parse_ascii_preview', {
                method: 'POST', headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({ filepath: filepath })
            }).then(r => r.json()).then(data => {
                if(data.error) { alert(data.error); return; }
                document.getElementById('ascii_preview').innerText = data.preview_text;
                
                ['col_ref', 'col_anc', 'col_glt', 'col_topo', 'col_brdf'].forEach(id => {
                    const sel = document.getElementById(id);
                    sel.innerHTML = '';
                    data.headers.forEach((h, i) => {
                        const opt = document.createElement('option');
                        opt.value = i; opt.text = `${h} (Col ${i+1})`;
                        sel.appendChild(opt);
                    });
                });
                
                // Auto-map naive logic based on names
                data.headers.forEach((h, idx) => {
                    let l = h.toLowerCase();
                    if(l.includes('ref') || l.includes('image')) document.getElementById('col_ref').value = idx;
                    if(l.includes('anc') || l.includes('obs')) document.getElementById('col_anc').value = idx;
                    if(l.includes('glt')) document.getElementById('col_glt').value = idx;
                    if(l.includes('topo')) document.getElementById('col_topo').value = idx;
                    if(l.includes('brdf')) document.getElementById('col_brdf').value = idx;
                });
            });
        }

        function applyAscii() {
            const filepath = document.getElementById('ascii_file').value;
            if(!filepath) return alert('Select an ASCII file first');
            
            const mapping = {
                ref: document.getElementById('col_ref').value,
                anc: document.getElementById('col_anc').value,
                glt: document.getElementById('col_glt').value,
                topo: document.getElementById('col_topo').value,
                brdf: document.getElementById('col_brdf').value
            };

            fetch('/api/apply_ascii', {
                method: 'POST', headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({ filepath: filepath, mapping: mapping })
            }).then(r => r.json()).then(data => {
                if(data.error) { alert(data.error); return; }
                listData['primary_files'] = data.ref;
                listData['anc_files'] = data.anc;
                listData['glt_files'] = data.glt;
                listData['topo_files'] = data.topo;
                listData['brdf_files'] = data.brdf;
                
                ['primary_files', 'anc_files', 'glt_files', 'topo_files', 'brdf_files'].forEach(updateListUI);
                document.getElementById('ascii_status').innerText = `✓ Auto-loaded ${data.ref.length} Reflectance files.`;
            });
        }

        // --- Payload Collection & API calls ---
        function collectPayload() {
            let chunk_r = document.getElementById('chunk_row').value;
            let chunk_c = document.getElementById('chunk_col').value;
            chunk_r = isNaN(chunk_r) ? chunk_r : parseInt(chunk_r);
            chunk_c = isNaN(chunk_c) ? chunk_c : parseInt(chunk_c);

            let corrections = [];
            if(document.getElementById('corr_topo').checked) corrections.push('topo');
            if(document.getElementById('corr_brdf').checked) corrections.push('brdf');

            return {
                timestamp: new Date().toISOString() + "Z",
                export_directory: document.getElementById('export_dir').value,
                inputs: {
                    images: listData['primary_files'],
                    enable_ascii_loader: document.getElementById('enable_ascii').checked,
                    anc_files: listData['anc_files'],
                    trait_files: listData['trait_files'],
                    enable_precomputed_correction: document.getElementById('enable_precomputed').checked,
                    topo_json_files: document.getElementById('enable_precomputed').checked ? listData['topo_files'] : [],
                    brdf_json_files: document.getElementById('enable_precomputed').checked ? listData['brdf_files'] : [],
                    use_glt: document.getElementById('enable_glt').checked,
                    glt_files: document.getElementById('enable_glt').checked ? listData['glt_files'] : []
                },
                parameters: {
                    data_type: document.getElementById('data_type').value,
                    data_type_export: document.getElementById('data_type_export').value,
                    bad_bands: document.getElementById('bad_bands').value,
                    corrections: corrections,
                    chunk_size: [chunk_r, chunk_c],
                    resample_method: document.getElementById('resample_method').value
                }
            };
        }

        function generatePreview() {
            fetch('/api/generate', {
                method: 'POST', headers: {'Content-Type': 'application/json'},
                body: JSON.stringify(collectPayload())
            }).then(r => r.json()).then(data => {
                let txt = "// Processed & Reordered JSON Config:\\n\\n";
                if(data.feedback && data.feedback.length > 0) {
                    txt += data.feedback.join("\\n\\n") + "\\n\\n";
                }
                txt += JSON.stringify(data.config, null, 4);
                document.getElementById('json_preview').innerText = txt;
            });
        }

        function saveJson() {
            const payload = collectPayload();
            const saveDir = document.getElementById('save_dir').value;
            const filename = document.getElementById('save_filename').value;
            
            if(!saveDir) return alert('Please select a JSON Save Folder first.');

            fetch('/api/save', {
                method: 'POST', headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({ payload: payload, dir: saveDir, filename: filename })
            }).then(r => r.json()).then(data => {
                let txt = "";
                if(data.feedback && data.feedback.length > 0) {
                    txt += data.feedback.join("\\n\\n") + "\\n\\n";
                }
                if(data.error) {
                    txt += "// ERROR: " + data.error + "\\n\\n";
                    alert('Error saving config.');
                } else {
                    txt += `// SUCCESS: File written to disk!\\n// Destination: ${data.path}\\n\\n`;
                    txt += JSON.stringify(data.config, null, 4);
                    alert(`Configuration successfully saved to:\\n${data.path}`);
                }
                document.getElementById('json_preview').innerText = txt;
            });
        }
        
        // Initialize default save dir visually
        window.onload = function() {
            fetch('/api/fs?path=.')
                .then(r => r.json())
                .then(data => document.getElementById('save_dir').value = data.path);
        };
    </script>
</body>
</html>
"""

@app.route('/')
def index():
    return HTML_TEMPLATE

@app.route('/api/fs', methods=['GET'])
def fs_browse():
    req_path = request.args.get('path', os.path.expanduser('~'))
    if req_path == '.':
        req_path = os.getcwd()

    target_path = os.path.abspath(req_path)

    if not os.path.exists(target_path):
        target_path = os.path.abspath(os.sep) # fallback to root

    dirs = []
    files = []

    try:
        for item in os.listdir(target_path):
            full_path = os.path.join(target_path, item)
            # Skip hidden files
            if item.startswith('.'):
                continue
            if os.path.isdir(full_path):
                dirs.append({'name': item, 'path': full_path})
            else:
                files.append({'name': item, 'path': full_path})
    except Exception as e:
        return jsonify({"error": str(e)}), 403

    dirs.sort(key=lambda x: x['name'].lower())
    files.sort(key=lambda x: x['name'].lower())

    return jsonify({
        'path': target_path,
        'parent': os.path.dirname(target_path),
        'dirs': dirs,
        'files': files
    })

@app.route('/api/parse_ascii_preview', methods=['POST'])
def parse_ascii_preview():
    filepath = request.json.get('filepath')
    if not filepath or not os.path.exists(filepath):
        return jsonify({"error": "Invalid file path"})

    try:
        with open(filepath, 'r', encoding='utf-8', errors='strict') as f:
            lines = [f.readline() for _ in range(50)]
        raw_preview = [l.strip('\r\n') for l in lines[:3] if l.strip()]

        header_line = raw_preview[0] if raw_preview else ""
        delimiter = ',' if ',' in header_line else '\t' if '\t' in header_line else ';' if ';' in header_line else None

        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=delimiter) if delimiter else csv.reader(f, delimiter=' ', skipinitialspace=True)
            all_rows = [row for row in reader if row and not row[0].startswith('#')]

        headers = [h.strip() for h in all_rows[0]] if all_rows else []
        return jsonify({
            "headers": headers,
            "preview_text": "\n".join(raw_preview)
        })
    except Exception as e:
        return jsonify({"error": str(e)})

@app.route('/api/apply_ascii', methods=['POST'])
def apply_ascii():
    data = request.json
    filepath = data.get('filepath')
    mapping = {k: int(v) for k, v in data.get('mapping', {}).items()}

    try:
        with open(filepath, 'r', encoding='utf-8', errors='strict') as f:
            header_line = f.readline()
            f.seek(0)
            delimiter = ',' if ',' in header_line else '\t' if '\t' in header_line else ';' if ';' in header_line else None
            reader = csv.reader(f, delimiter=delimiter) if delimiter else csv.reader(f, delimiter=' ', skipinitialspace=True)
            all_rows = [row for row in reader if row and not row[0].startswith('#')]

        data_rows = [[cell.strip() for cell in row] for row in all_rows[1:]]

        results = {'ref': [], 'anc': [], 'glt': [], 'topo': [], 'brdf': []}

        for row in data_rows:
            if mapping.get('ref') < len(row) and row[mapping['ref']]: results['ref'].append(row[mapping['ref']])
            if mapping.get('anc') < len(row) and row[mapping['anc']]: results['anc'].append(row[mapping['anc']])
            if mapping.get('glt') < len(row) and row[mapping['glt']]: results['glt'].append(row[mapping['glt']])
            if mapping.get('topo') < len(row) and row[mapping['topo']]: results['topo'].append(row[mapping['topo']])
            if mapping.get('brdf') < len(row) and row[mapping['brdf']]: results['brdf'].append(row[mapping['brdf']])

        return jsonify(results)
    except Exception as e:
        return jsonify({"error": str(e)})

@app.route('/api/generate', methods=['POST'])
def generate():
    raw_payload = request.json
    processed_config, feedback_message = process_and_reorder_config(raw_payload)
    return jsonify({
        "config": processed_config,
        "feedback": feedback_message
    })

@app.route('/api/save', methods=['POST'])
def save():
    data = request.json
    raw_payload = data.get('payload')
    save_dir = data.get('dir')
    filename = data.get('filename', 'local_job_config.json').strip()

    if not filename.endswith('.json'):
        filename += '.json'

    processed_config, feedback_message = process_and_reorder_config(raw_payload)

    if len(feedback_message) > 0:
        return jsonify({"error": "File is not saved due to warnings.", "feedback": feedback_message})

    full_path = os.path.join(save_dir, filename)
    try:
        with open(full_path, 'w') as f:
            json.dump(processed_config, f, indent=4)
        return jsonify({
            "path": full_path,
            "config": processed_config,
            "feedback": feedback_message
        })
    except Exception as e:
        return jsonify({"error": str(e)})

if __name__ == '__main__':
    if len(sys.argv)>1:
        port_user = int(sys.argv[1])
    elif len(sys.argv)==1:
        port_user = 5005    
    #app.run(host='127.0.0.1', port=port_user, debug=True)
    app.run(host='0.0.0.0', port=port_user, debug=True)

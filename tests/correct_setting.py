
BAD_BAND_RANGE =  [[300,415],[900,1430],[1800,2600]]  # [[300,405],[1337,1430],[1800,1960],[2450,2600]]

glint_setting_gao_1 = {
      "type": "gao",
      "correction_wave": 1640,
      "apply_mask": [
         [
            "ndi",
            {
               "band_1": 550,
               "band_2": 2150,
               "min": 0,
               "max": 1
            }
         ],
         [
            "ndi",
            {
               "band_1": 850,
               "band_2": 660,
               "min": -1,
               "max": 0.1
            }
         ]
      ],
      "truncate":False

}

glint_setting_hochberg = {
      "type": "hochberg",
      "correction_wave": 2150,
      "apply_mask": [
         [
            "ndi",
            {
               "band_1": 550,
               "band_2": 2150,
               "min": 0,
               "max": 1
            }
         ]
      ],
      "truncate":False

}


glint_setting_gao = {
      "type": "gao",
      "correction_wave": 1640,
      "apply_mask": [
         [
            "ancillary",
            {
               "name": "slope",
               "min": "-inf",
               "max": 0.001
            }
         ],
         [
            "ndi",
            {
               "band_1": 550,
               "band_2": 2150,
               "min": 0,
               "max": 1
            }
         ],
      ],
      "truncate":False
}

topo_dict = {
      "type": "scs+c",
      "calc_mask": [
         [
            "ndi",
            {
               "band_1": 850,
               "band_2": 660,
               "min": 0.05,
               "max": 1.0
            }
         ],
         [
            "ancillary",
            {
               "name": "slope",
               "min": 0.08726646259971647,
               "max": "+inf"
            }
         ],
         [
            "ancillary",
            {
               "name": "cosine_i",
               "min": 0.12,
               "max": "+inf"
            }
         ]
      ],
      "apply_mask": [
         [
            "ndi",
            {
               "band_1": 850,
               "band_2": 660,
               "min": 0.05,
               "max": 1.0
            }
         ],
         [
            "ancillary",
            {
               "name": "slope",
               "min": 0.08726646259971647,
               "max": "+inf"
            }
         ],
         [
            "ancillary",
            {
               "name": "cosine_i",
               "min": 0.12,
               "max": "+inf"
            }
         ]
      ],
      "c_fit_type": "nnls",
}

brdf_dict = {
      "solar_zn_type": "line",
      "type": "flex",
      "geometric": "li_dense_r",
      "volume": "ross_thick",
      "b/r": 2.5,
      "h/b": 2,
      "sample_perc": 0.5,
      "interp_kind": "linear",
      "calc_mask": [
         [
            "ndi",
            {
               "band_1": 850,
               "band_2": 660,
               "min": 0.05,
               "max": 1.0
            }
         ],
         [
            "kernel_finite",
            {}
         ],
         [
            "ancillary",
            {
               "name": "slope",
               "min": 0.01,
               "max": "+inf"
            }
         ],
         [
            "ancillary",
            {
               "name": "sensor_zn",
               "min": 0.03490658503988659,
               "max": "inf"
            }
         ],
      ],
      "apply_mask": [
         [
            "ancillary",
            {
               "name": "slope",
               "min": 0.01,
               "max": "+inf"
            }
         ],
         [
            "ndi",
            {
               "band_1": 850,
               "band_2": 660,
               "min": 0.05,
               "max": 1.0
            }
         ]
      ],
      "bin_type": "dynamic",
      "num_bins": 18,
      "ndvi_bin_min": 0.05,
      "ndvi_bin_max": 1.0,
      "ndvi_perc_min": 10,
      "ndvi_perc_max": 95
}
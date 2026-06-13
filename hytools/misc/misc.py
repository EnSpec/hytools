# -*- coding: utf-8 -*-
"""
HyTools:  Hyperspectral image processing library
Copyright (C) 2021 University of Wisconsin

Authors: Adam Chlus, Zhiwei Ye, Philip Townsend.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
from itertools import tee
import numpy as np
import pandas as pd


def progbar(curr, total, full_progbar = 100):
    '''Display progress bar.

    Gist from:

    https://gist.github.com/marzukr/3ca9e0a1b5881597ce0bcb7fb0adc549

    Args:
        curr (int, float): Current task level.
        total (int, float): Task level at completion.
        full_progbar (TYPE): Defaults to 100.

    Returns:
        None.

    '''
    frac = curr/total
    filled_progbar = round(frac*full_progbar)
    print('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')


def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def set_brdf(hy_obj,brdf_dict):
    hy_obj.brdf = brdf_dict

def set_topo(hy_obj,topo_dict):
    hy_obj.topo = topo_dict

def update_brdf(hy_obj,args):
    hy_obj.brdf[args['key']] = args['value']

def update_topo(hy_obj,args):
    hy_obj.topo[args['key']] = args['value']

def set_glint(hy_obj,glint_dict):

    # If the type is hedley, need to specify deep water area
    if glint_dict['type'] == 'Hedley':
        glint_dict['deep_water_sample'] = glint_dict['deep_water_sample'][hy_obj.file_name]

    hy_obj.glint = glint_dict

def update_topo_group(subgroup_dict_in):

    subgroup_dict = {}
    group_tag_list=[]

    for file_name in subgroup_dict_in.keys():
        group_tag = subgroup_dict_in[file_name]
        if group_tag in subgroup_dict:
            subgroup_dict[group_tag]+=[file_name]
        else:
            subgroup_dict[group_tag]=[file_name]
            group_tag_list+=[group_tag]

    update_name_list=[]
    for group_tag in subgroup_dict.keys():
        update_name_list+=[subgroup_dict[group_tag]]

    return update_name_list,group_tag_list

# subsample non-zero pixels in a class mask using regular grid sampling 
# and keep sampling locations to a CSV file
def regular_grid_sampling_class_mask(class_mask, sample_perc, csv_path = None):

    # safety check
    if not(0 < sample_perc <=1):
        raise ValueError("sample_perc must be in the range of (0, 1]. ")
    
    # create the np.ndarray to store the output mask
    class_mask_out = class_mask.copy()

    #  calcalate the number of valid pixels
    valid_mask = class_mask_out != 0
    n_valid = np.sum(valid_mask)

    if n_valid == 0: 
        if csv_path is not None: 
            pd.DataFrame(columns = ["rowID", "colID", "class_value"]).to_csv(
                csv_path, index = False
            )
            return class_mask_out.astype(np.int8)
        
    if sample_perc == 1:
        keep_mask = valid_mask
        grid_step = 1

    else: 
        # calculate the number of pixels to sample
        n_target = int(n_valid * sample_perc)

        row_grid, col_grid = np.indices(class_mask_out.shape)

        max_grid_step = max(class_mask_out.shape)

        best_keep_mask = None
        best_grid_step = None
        best_diff = np.inf
        keep_n_step = None

        for candidate_step in range(1, max_grid_step + 1):
            
                keep_mask_candidate = (
                    valid_mask
                    & ((row_grid) % grid_step == 0)
                    & ((col_grid) % grid_step == 0)
                )

                n_keep = np.sum(keep_mask_candidate)
                diff = abs(n_keep - n_target)
                
                # print out the grid_step, n_keep, and diff number
                

                if diff < best_diff:
                    best_diff = diff
                    best_keep_mask = keep_mask_candidate
                    best_grid_step = candidate_step
                    keep_n_step = n_keep

                # Exact match, no need to keep seraching
                if diff == 0: 
                    break

        keep_mask = best_keep_mask
        grid_step = best_grid_step

    if csv_path is not None: 
        keep_rows, keep_cols = np.where(keep_mask)

        sampling_df = pd.DataFrame({
            "rowID": keep_rows,
            "colID": keep_cols,
            "class_value": class_mask_out[keep_rows, keep_cols],
            "grid_Step": grid_step
        })

        sampling_df.to_csv(csv_path, index=False)

    # this is equivalent to sampling (1-sample_perc) pixels and setting them to 0
    class_mask_out[valid_mask & (~keep_mask)] = 0

    return class_mask_out.astype(np.int8)







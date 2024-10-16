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

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
import re

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
    """Re-organize subgroup dictionary.
        Input dictionary use group name as key,
        and the return is two lists of paired image names and group names.  

    Args:
        subgroup_dict_in (dictionary): subgroup info 

    Returns:
        tuple: A tuple containing:
            - update_name_list (list): name_list.
            - group_tag_list (list): group_tag_list.

    """

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

def wkt_parse_update(wkt_string):
    '''Parse and update wkt string for map projection, in case "unnamed" projection is provided.

    Args:
        wkt_string (string): wkt string.

    Returns:
        new_wkt_string (string): updated wkt string.
        zone_tag (string): zone id and N/S
    '''

    # 1. Extract Projection / CRS Name (First quoted string in PROJCS, GEOGCS
    proj_name_match = re.search(
        r'(?:PROJCS|GEOGCS|COMPOUNDCRS|PROJCRS)\s*\[\s*"([^"]+)"', wkt_string
    )
    projection_name = proj_name_match.group(1) if proj_name_match else None

    # 2. Extract PROJECTION and value pair
    proj_type = re.findall(
        r'PROJECTION\s*\[\s*"([^"]+)"\s*\]', wkt_string
    )

    # 3. Extract all PARAMETER name and value pairs
    parameters = re.findall(
        r'PARAMETER\s*\[\s*"([^"]+)"\s*,\s*([-\d.]+)\s*\]', wkt_string
    )
    param_dict = {param[0]: float(param[1]) for param in parameters}

    # 4. Extract EPSG Code (Looks for AUTHORITY["EPSG", "XXXX"] or ID["EPSG", "XXXX"])
    epsg_matches = re.findall(
        r'(?:AUTHORITY|ID)\s*\[\s*"EPSG"\s*,\s*"(\d+)"\s*\]',
        wkt_string,
        re.IGNORECASE,
    )
    # The main EPSG code of the CRS is typically the last AUTHORITY block in a WKT string
    main_epsg = epsg_matches[-1] if epsg_matches else None

    if "Transverse_Mercator" == proj_type[0]:
        if "utm zone" in projection_name.lower():

            return wkt_string, projection_name.split('UTM zone ')[1]  # do not update

        # else: # 'unnamed' or something else
        zone_id = int((param_dict["central_meridian"]-3 +180)/6)+1
        new_wkt_string = wkt_string
        if param_dict["false_northing"]==1e7:
            zone_tag = str(zone_id) + 'S'
            main_epsg = 32700 + zone_id
        elif param_dict["false_northing"]==0.0:
            zone_tag = str(zone_id) + 'N'
            main_epsg = 32600 + zone_id
        else:
            return wkt_string, None # do not update

        # Rename projection_name and append EPSG code
        new_wkt_string = new_wkt_string.replace(projection_name, f"WGS 84 / UTM zone {zone_id}", 1)
        new_wkt_string = f'AXIS["Northing",NORTH],AUTHORITY["EPSG","{main_epsg}"]'.join(new_wkt_string.rsplit('AXIS["Northing",NORTH]', 1))
        return new_wkt_string, zone_tag
    #else:
    return wkt_string, None  # do not update

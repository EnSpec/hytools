# -*- coding: utf-8 -*-
"""
This module contains functions to apply the Modified topographic correction (SCS+C)
described in the following paper:

Richter, R., Kellenberger, T., & Kaufmann, H. (2009).
Comparison of topographic correction methods.
Remote Sensing, 1(3), 184-196.
https://doi.org/10.3390/rs1030184

Topographic correction consists of the following steps:


"""
import numpy as np

def calc_modminn_coeffs(hy_obj,topo_dict):
    '''

    Args:
        hy_obj (TYPE): DESCRIPTION.

    Returns:
        None.

    '''
    hy_obj.topo =topo_dict

    cos_i = hy_obj.cosine_i()
    i = np.rad2deg(np.arccos(cos_i))
    solar_zn = hy_obj.get_anc('solar_zn',radians=False)

    solar_zn_t = np.zeros(solar_zn.shape)
    solar_zn_t[solar_zn < 45] = solar_zn[solar_zn < 45] +20
    solar_zn_t[(solar_zn >= 45) & (solar_zn <= 55)] = solar_zn[(solar_zn >= 45) & (solar_zn <= 55)] +15
    solar_zn_t[solar_zn > 55] = solar_zn[solar_zn > 55] +10

    #Create NDVI mask to seperate vegetation
    ir = hy_obj.get_wave(850)
    red = hy_obj.get_wave(660)
    ndvi = (ir-red)/(ir+red)
    veg_mask = ndvi > 0.2

    c_factors = np.ones((2,hy_obj.lines,hy_obj.columns))
    c_factors[:] = cos_i/np.cos(np.radians(solar_zn_t))

    # Non vegetation correction factor
    c_factors[0][~veg_mask] =  c_factors[0][~veg_mask]**(1/2)
    c_factors[1][~veg_mask] =  c_factors[1][~veg_mask]**(1/2)

    # Vegetation correction factor
    c_factors[0][veg_mask] =  c_factors[0][veg_mask]**(3/4)
    c_factors[1][veg_mask] =  c_factors[1][veg_mask]**(1/3)

    #Adjust correction factors to prevent too strong correction
    c_factors[c_factors <.25] = .25
    c_factors[c_factors > 1] = 1

    #Correct pixels only where i > threshold
    c_factors[0][i < solar_zn_t] = 1
    c_factors[1][i < solar_zn_t] = 1

    c_factors[0][ir == hy_obj.no_data] = 1
    c_factors[1][ir == hy_obj.no_data] = 1

    hy_obj.ancillary['mm_c_factor'] = c_factors

def apply_modminn(hy_obj,data,dimension,index):
    ''' Apply SCSS correction to a slice of the data

    Args:
        hy_obj (TYPE): DESCRIPTION.
        band (TYPE): DESCRIPTION.
        index (TYPE): DESCRIPTION.

    Returns:
        band (TYPE): DESCRIPTION.

    '''

    if 'mm_c_factor' not in hy_obj.ancillary.keys():
        calc_modminn_coeffs(hy_obj)

    #Convert to float
    data = data.astype(np.float32)

    wave_mask =hy_obj.wavelengths >=720

    if dimension == 'line':
        #index= 3000
        #data = hy_obj.get_line(3000)
        data[:,wave_mask] = data[:,wave_mask]*hy_obj.ancillary['mm_c_factor'][1,index,:][:,np.newaxis]
        data[:,~wave_mask] = data[:,~wave_mask]*hy_obj.ancillary['mm_c_factor'][0,index,:][:,np.newaxis]

    elif dimension == 'column':
        #index= 300
        #data = hy_obj.get_column(index)
        data[:,wave_mask] = data[:,wave_mask]*hy_obj.ancillary['mm_c_factor'][1,:,index][:,np.newaxis]
        data[:,~wave_mask] = data[:,~wave_mask]*hy_obj.ancillary['mm_c_factor'][0,:,index][:,np.newaxis]

    elif dimension == 'band':
        #index= 50
        #data = hy_obj.get_band(index)
        if hy_obj.wavelengths[index] >=720:
            cf_index = 1
        else:
            cf_index = 0
        data = data * hy_obj.ancillary['mm_c_factor'][cf_index]

    elif dimension == 'chunk':
        #index = 200,501,3000,3501
        x1,x2,y1,y2 = index
        #data = hy_obj.get_chunk(x1,x2,y1,y2)
        data[:,:,wave_mask]  = data[:,:,wave_mask]*hy_obj.ancillary['mm_c_factor'][1,y1:y2,x1:x2][:,:,np.newaxis]
        data[:,:,~wave_mask] = data[:,:,~wave_mask]*hy_obj.ancillary['mm_c_factor'][0,y1:y2,x1:x2][:,:,np.newaxis]

    elif dimension == 'pixels':
        #index = [[2000,2001],[200,501]]
        y,x = index
        #data = hy_obj.get_pixels(y,x)
        data[:,wave_mask] = data[:,wave_mask]*hy_obj.ancillary['mm_c_factor'][1,y,x][:, np.newaxis]
        data[:,~wave_mask] = data[:,~wave_mask]*hy_obj.ancillary['mm_c_factor'][0,y,x][:, np.newaxis]
    return data

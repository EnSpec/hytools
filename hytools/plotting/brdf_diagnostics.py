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

Plotting functions for BRDF
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def universal_diagno_plot(hy_obj,config_dict):
    ''' Generate a diagnostic plot of BRDF correction results.
    '''
    #Flip sign of zenith angle at minimum
    sensor_zn = hy_obj.get_anc('sensor_zn',radians =False)
    sensor_zn[~hy_obj.mask['no_data']] = np.nan
    for i,line in enumerate(sensor_zn):
        line[:np.nanargmin(line)] *= -1
        sensor_zn[i] = line
    sensor_zn = (sensor_zn[hy_obj.mask['calc_brdf']]//2)*2

    diagno_df = pd.DataFrame()
    diagno_df['sensor_zn'] =sensor_zn

    bands = [hy_obj.wave_to_band(wave) for wave in config_dict['brdf']['diagnostic_waves']]
    for band_num in bands:
        band = hy_obj.get_band(band_num,mask='calc_brdf')
        diagno_df['uncorr_%s' % band_num] =  band

        band = hy_obj.get_band(band_num,
                               corrections = hy_obj.corrections + ['brdf'],
                               mask='calc_brdf')

        diagno_df['corr_%s' % band_num] =  band
        fvol, fgeo, fiso  = hy_obj.brdf['coeffs'][band_num]

        brdf = fvol*hy_obj.ancillary['k_vol']
        brdf += fgeo*hy_obj.ancillary['k_geom']
        brdf+=fiso
        brdf = brdf[hy_obj.mask['calc_brdf']]
        diagno_df['brdf_%s' % band_num] =  brdf

    # Average every 2 degrees of zenith angle
    diagno_df=  diagno_df.groupby(by= 'sensor_zn').mean()

    fig = plt.figure(figsize= (8,6))
    fig.suptitle(hy_obj.base_name)
    for a,band_num in enumerate(bands,start=1):
        ax = fig.add_subplot(2,2,a)
        ax.plot(diagno_df.index,diagno_df['brdf_%s' % band_num],c='k',ls ='--')
        ax.scatter(diagno_df.index,diagno_df['uncorr_%s' % band_num],marker ='o',fc='w',ec='k')
        ax.scatter(diagno_df.index,diagno_df['corr_%s' % band_num],marker ='o',fc='k',ec='k')
        ax.text(.85,.9, "%s nm" % int(hy_obj.wavelengths[band_num]), transform=ax.transAxes,
                ha = 'center', fontsize = 12)
        if a > 2:
            ax.set_xlabel('View zenith angle')
        if a in [1,3]:
            ax.set_ylabel('Reflectance')

    #Create legend
    custom_points = []
    custom_points.append(Line2D([0],[0], marker = 'o',label='Uncorrected',
                          markerfacecolor='w', markersize=10,lw=0,markeredgecolor='k'))
    custom_points.append(Line2D([0],[0], marker = 'o',label='Corrected',
                          markerfacecolor='k', markersize=10,lw=0,markeredgecolor='k'))
    custom_points.append(Line2D([0],[1],label='Modeled BRDF',c='k', ls ='--'))
    ax.legend(handles=custom_points, loc='center',frameon=False,
              bbox_to_anchor=(-.15, -.3), ncol =3,columnspacing = 1.5,labelspacing=.25)

    plt.savefig("%s%s_brdf_plot.png" % (config_dict['export']['output_dir'],hy_obj.base_name),
                bbox_inches = 'tight')
    plt.close()

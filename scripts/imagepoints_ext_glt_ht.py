
import os, sys
import numpy as np
import argparse
import pandas as pd

try:
    from osgeo import gdal, osr, ogr
    has_gdal=True
except ModuleNotFoundError:
    print("No gdal installed")
    has_gdal=False

import hytools as ht
from hytools.misc.point import local_point2spec, subset_band_list #*
#warnings.filterwarnings("ignore")
#np.seterr(divide='ignore', invalid='ignore')


def obs_point2spec(hyObj, img_row, img_col):

    b1=hyObj.get_anc('sensor_az')[ img_row,img_col]
    b2=hyObj.get_anc('sensor_zn')[ img_row,img_col]
    b3=hyObj.get_anc('solar_az')[ img_row,img_col]
    b4=hyObj.get_anc('solar_zn')[ img_row,img_col]
    b5=hyObj.get_anc('slope')[ img_row,img_col]
    b6=hyObj.get_anc('aspect')[img_row,img_col]

    obs_data = np.vstack((img_row,img_col,b1,b2,b3,b4,b5,b6))

    obs_df = pd.DataFrame(obs_data.T, columns=['img_row','img_col','sensor_az','sensor_zn','solar_az','solar_zn','slope','azimuth'])
    return obs_df

def rasterize_polygon(hyObj,polygon_fn, key_id,use_glt_bool):
    ''' Rasterize polygon based on image georeference and boundary
    
    '''
    source_ds = ogr.Open(polygon_fn)

    source_layer = source_ds.GetLayer()

    field_list = []
    ldefn = source_layer.GetLayerDefn()
    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        field_list.append(fdefn.name)
    print(field_list[1])

    if not (key_id in field_list):
        print('Field "',key_id,'" is not in the shapefile!')
        return (None, None)

    tmp_mem_driver=ogr.GetDriverByName('MEMORY')

    dest = tmp_mem_driver.CreateDataSource('tempData')

    mem_lyr = dest.CopyLayer(source_layer,'newlayer',['OVERWRITE=YES'])
    FeatureCount= mem_lyr.GetFeatureCount()
    # Add a new field
    new_field = ogr.FieldDefn('tempFID', ogr.OFTInteger)
    mem_lyr.CreateField(new_field)

    lookup_dict={}
    for i, feature in enumerate(mem_lyr):

        feature.SetField('tempFID', i+1)  # key step1
        lookup_dict[str(i+1)]=feature.GetField(key_id)
        mem_lyr.SetFeature(feature)  # key step 2

    if use_glt_bool:
        out_col = hyObj.columns_glt
        out_row = hyObj.lines_glt
        out_transform = hyObj.glt_transform
        out_proj = hyObj.glt_projection
    else:
        out_col = hyObj.columns
        out_row = hyObj.lines
        out_transform = hyObj.transform
        out_proj = hyObj.projection


    if FeatureCount<255:
        target_ds = gdal.GetDriverByName('MEM').Create('', out_col, out_row, 1, gdal.GDT_Byte)

        nodata_val=255
    else:
        if (FeatureCount>255 and FeatureCount < 32767):
            target_ds = gdal.GetDriverByName('MEM').Create('', out_col, out_row, 1, gdal.GDT_Int16)
            nodata_val=-9999
        else:    # >32767
            target_ds = gdal.GetDriverByName('MEM').Create('', out_col, out_row, 1, gdal.GDT_Int32)
            nodata_val=-9999

    target_ds.SetGeoTransform(out_transform)

    target_ds.SetProjection(out_proj)
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata_val)

    gdal.RasterizeLayer(target_ds, [1], mem_lyr, options=["ATTRIBUTE=tempFID" ,"ALL_TOUCHED=FALSE"])

    return (target_ds,lookup_dict)

def gen_df_polygon(hyObj, target_ds, lookup_dict, imgsrs2latlon, uid, use_glt_bool):
    '''Generate a dataframe that stores the location, UID of all the points within the polygons 
    
    Parameters
    ----------

    hyObj:            HyTools file object
    target_ds:      GDAL raster dataset
                            one band raster in which each polygon has unique digital number
    lookup_dict:    dictionary
                            a distionary linking polygon DN in raster (target_ds) and the UID in the polygons attribute table
    imgsrs2latlon: coordinate transformation object
                            transform from georeferenced coordinates of the image to LAT LON
    uid:                 str
                            the user specified unique polygon ID name from the attribute table of the shapefile

    Returns
    -------
    return_df:   pandas dataframe
                        a dataframe that stores the location, UID of all the points within the polygons
    
    '''

    poly_raster=target_ds.GetRasterBand(1).ReadAsArray()

    data_type=target_ds.GetRasterBand(1).DataType

    if data_type == gdal.GDT_Byte:
        ind=np.where((poly_raster>0) & (poly_raster<255)  )
    else:
        if data_type == gdal.GDT_Int16:
            ind=np.where((poly_raster>0) & (poly_raster<32767)  )
        else:
            ind=np.where(poly_raster>0  )

    total_point=len(ind[1])

    print(total_point,' points')

    if total_point==0:
        # polygons are not intersecting the image
        print( "No intersection.")
        return None

    return_df = pd.DataFrame(columns=['new_uid',uid,'img_col_glt','img_row_glt','img_col_raw','img_row_raw','lon','lat'])
    return_df = return_df.fillna(0) # with 0s rather than NaNs

    if use_glt_bool:
        ul_x, new_x_resolution, new_x_rot, ul_y, new_y_rot, new_y_resolution = hyObj.glt_transform
    else:
        ul_x, new_x_resolution, new_x_rot, ul_y, new_y_rot, new_y_resolution = hyObj.transform

    sub_id_dict = {}
    for key in lookup_dict.keys():
        sub_id_dict[key]=0

    # add polygon ID, and point order number within the same polygon
    for index in range(total_point):

        row=ind[0][index]
        col=ind[1][index]

        if use_glt_bool:
            row_post_glt = hyObj.glt_y[row,col] - 1
            col_post_glt = hyObj.glt_x[row,col] - 1 # zero-based
        else:
            row_post_glt = row
            col_post_glt = col

        poly_id=poly_raster[row,col]
        poly_id_code=lookup_dict[str(poly_id)]

        x_coord = ul_x + (col+0.5)*new_x_resolution + (row+0.5)*new_x_rot
        y_coord = ul_y + (col+0.5)*new_y_rot +            (row+0.5)*new_y_resolution

        lat, lon, _ = imgsrs2latlon.TransformPoint(x_coord, y_coord) #lon, lat, _ = imgsrs2latlon.TransformPoint(x_coord, y_coord)

        sub_id = sub_id_dict[str(poly_id)]
        sub_id_dict[str(poly_id)]+=1

        temp_df = pd.DataFrame([['{}_{}'.format(poly_id_code, sub_id),poly_id_code, col,row,col_post_glt,row_post_glt, lon,lat]], columns=['new_uid', uid, 'img_col_glt','img_row_glt','img_col_raw','img_row_raw','lon','lat'])

        return_df = pd.concat([return_df,temp_df],ignore_index=True)  #return_df.append(temp_df,ignore_index=True)


    return return_df

def local_polygon2spec(hyObj, poly_shp, uid, use_band_list=False, band_list=[],use_glt_bool=False):
    """Extract spectra with points within the boundary of polygons from the hyperspectral image
    
    Steps:
    1, Rasterize polygon based on image georeference
    2, Get locations of the points of interest from the raster
    3, Overlapping points and the hyperspectral image, and extract spectra  
    
    Parameters
    ----------
    hyObj :              HyTools file object
    poly_shp:          str
                            full filename of the polygon shapefile 
    uid:                   str
                            the user specified unique polygon ID name from the attribute table of the shapefile
    use_band_list: boolean
                            default True; whether to use a subset of bands
    band_list:         list or numpy array                         
                            default is a blank list
                            if it is a list, it should be one like [5,6,7,8,9, 12]
                            if it is a numpy array, it should be the same size as hyObj.bad_bands with only True or False in the array
    use_glt_bool: boolean
                        default False;                        
            
    Returns
    -------
    point_df: pandas dataframe
                    it include all the location and spectra information for all points within the polygons
    
    """

    if use_glt_bool:
        img_srs = osr.SpatialReference(wkt=hyObj.glt_projection)
    else:
        img_srs = osr.SpatialReference(wkt=hyObj.projection)

    latlon_wgs84 = osr.SpatialReference()
    latlon_wgs84.ImportFromEPSG ( 4326 )

    # LAT LON will be the only georeferenced coordinates kept in the result
    imgsrs2latlon = osr.CoordinateTransformation (img_srs, latlon_wgs84)

    # convert polygon geometry into raster with the same size of the image, and store UID in a lookup dictionary
    target_ds, lookup_dict=rasterize_polygon(hyObj,poly_shp,uid,use_glt_bool)

    if target_ds is None:
        return None

    # generate a dataframe that stores the location, UID of all the points within the polygons
    point_df = gen_df_polygon(hyObj, target_ds, lookup_dict, imgsrs2latlon, uid, use_glt_bool)

    if point_df is None :
        return None

    # extract full spectra information from image based on points locations
    #spec_data = extract_from_point(hyObj, point_df)

    spec_data = hyObj.get_pixels(point_df['img_row_raw'].values.astype(np.int16),point_df['img_col_raw'].values.astype(np.int16))

    # determine the column names of the spectra dataframe based on wavelengths
    if hyObj.wavelength_units.lower()[:4]=='micr':
        new_band_name = ['B{:0.3f}'.format(x) for x in hyObj.wavelengths]
    elif hyObj.wavelength_units.lower()[:4]=='nano' :
        new_band_name = ['B{:04d}'.format(int(x)) for x in hyObj.wavelengths]
    else:
        new_band_name = ['B{:d}'.format(x+1) for x in range(hyObj.bands)]

    if hyObj.file_type in ['ncav']:
        spec_df = pd.DataFrame(spec_data, columns=new_band_name) #spec_df = pd.DataFrame(spec_data.T, columns=new_band_name)
    else:
        spec_df = pd.DataFrame(spec_data, columns=new_band_name)

    # perform the subsetting of the columns in the dataframe according to the band_list or hyObj.bad_bands
    spec_df = subset_band_list(hyObj,spec_df,use_band_list, band_list)

    # merge location information and spectra information
    point_df = pd.concat([point_df,spec_df], axis=1, join='inner')

    return point_df

def main():
    parser = argparse.ArgumentParser(description='Export fractional cover image by EndMember csv')

    parser.add_argument('-i', type=str, required=True,help='Input image file name')
    parser.add_argument('-pnt', type=str, required=True,help='CSV filename or shapefile')
    parser.add_argument('-od', type=str, required=True,help='Output folder')
    parser.add_argument('-uid', type=str, required=True,help='Unique ID in the vector file')
    parser.add_argument('-epsg', type=str, required=False, help='UTM EPSG code')

    parser.add_argument('-anc', type=str, required=False, help='Ancillary file / OBS file')
    parser.add_argument('-glt', type=str, default=None, required=False, help='External GLT ENVI file')

    parser.add_argument('-dt', type=str, default='envi', required=False, help="Data type of the image (default 'envi') ['envi','emit','ncav']", choices=['envi','emit','ncav'])

    parser.add_argument('-nnb', type=int, required=False, default=4,help='How many neighbors in the image should be sampled from the center', choices=[0,4,8])



    args = parser.parse_args()

    in_image_file = args.i
    out_path = args.od
    pnt_file = args.pnt
    uid = args.uid
    epsg_code = args.epsg
    n_neighbor_chose = args.nnb

    file_format = args.dt

    if not args.glt is None:
        glt_dict = {
            "glt_x": [args.glt,1],
            "glt_y": [args.glt,0]
        }
    else:
        glt_dict = {}

    if args.anc: #not args.anc is None
        anc_dict = {
         "path_length": [
            args.anc,
            0
         ],
         "sensor_az": [
            args.anc,
            1
         ],
         "sensor_zn": [
            args.anc,
            2
         ],
         "solar_az": [
            args.anc,
            3
         ],
         "solar_zn": [
            args.anc,
            4
         ],
         "phase": [
            args.anc,
            5
         ],
         "slope": [
            args.anc,
            6
         ],
         "aspect": [
            args.anc,
            7
         ],
         "cosine_i": [
            args.anc,
            8
         ],
         "utc_time": [
            args.anc,
            9
         ]
        }
    else:
        anc_dict = None

    if pnt_file.endswith('csv'):
        if args.epsg is None:
            epsg_code = None
        else:
            epsg_code = args.epsg


    hy_obj = ht.HyTools()
    hy_obj.read_file(in_image_file,file_format, glt_path=glt_dict, anc_path =anc_dict)

    lookup_glt_bool = False
    if file_format=='emit':
        lookup_glt_bool=True
    else:
        if glt_dict:
            lookup_glt_bool=True

    if has_gdal and (pnt_file.endswith('.shp') or pnt_file.endswith('.geojson') or pnt_file.endswith('.json')):
        out_df = local_polygon2spec(hy_obj, pnt_file, uid, use_band_list=False, band_list=[],use_glt_bool=lookup_glt_bool)
    else:
        if not pnt_file.endswith('.csv'):
            print("Point location file is not in CSV format")
            return

        pnt_df = pd.read_csv(pnt_file)
        if 'x_coord' in pnt_df.columns and 'y_coord' in pnt_df.columns:
            out_df = local_point2spec(hy_obj, pnt_file, uid, 'x_coord', 'y_coord', epsg_code, n_neighbor=n_neighbor_chose, use_band_list=False, band_list=[],use_glt_bool=lookup_glt_bool)
        elif 'lat' in pnt_df.columns and 'lon' in pnt_df.columns:
            print('latlon inside point file.')
            out_df = local_point2spec(hy_obj, pnt_file, uid, 'lon', 'lat', epsg_code, n_neighbor=n_neighbor_chose, use_band_list=False, band_list=[],use_glt_bool=lookup_glt_bool)
        else:
            print("Unknown coordinates column names")
            return

    if out_df is None:
        return

    img_base_name=os.path.basename(in_image_file).split('.')[0]

    out_df.insert(loc = 1,
          column = 'flightline',
          value = [img_base_name.split('_')[0]]*out_df.shape[0])

    out_df.to_csv(out_path+img_base_name+"_spec_df_asvc.csv",index=False)

    if not args.anc is None:
        img_row = out_df['img_row_raw'].values.astype(np.int16)
        img_col = out_df['img_col_raw'].values.astype(np.int16)
        out_obs_df = obs_point2spec(hy_obj, img_row, img_col)
        out_obs_df.insert(loc = 0,
          column = 'flightline',
          value = [img_base_name.split('_')[0]]*out_df.shape[0])
        out_obs_df.insert(loc = 0,
          column = 'new_uid',
          value = out_df['new_uid'])
        out_obs_df.to_csv(out_path+img_base_name+"_obs_df.csv",index=False)

if __name__== "__main__":
    main()



import pandas as pd
from .geog_utm import *

def local_transform_all_point(mapobj, point_df, uid, xcoord, ycoord,point_epsg_code):
    ''' Create a dataframe with image georeferenced coordinates of all points of interest
  
    '''

    if point_epsg_code is None:
        print("Default latlon")
        re_df = pd.DataFrame(point_df[[uid,xcoord,ycoord]])
        re_df.columns = [uid,'img_x','img_y']
        return re_df
    else:
        ycoord_arr = point_df[ycoord]
        xcoord_arr = point_df[xcoord]
        lat_arr,lon_arr=mapobj.convert_latlon(xcoord_arr,ycoord_arr)
        re_df = point_df[[uid, xcoord, ycoord]].join(pd.DataFrame(np.array((lat_arr,lon_arr)).T))
        re_df.columns=[uid,'img_x','img_y','lat','lon']
        return re_df

def get_neighbor(hyObj, point_coord_df, n_neighbor, uid, point_epsg_code,mapobj,use_glt_bool):
    ''' Create a dataframe with columns and lines of all image space neighbors of points of interest
    
    '''
    if use_glt_bool:
        ul_x, new_x_resolution, new_x_rot, ul_y, new_y_rot, new_y_resolution = hyObj.glt_transform
        print(hyObj.glt_projection,hyObj.glt_map_info)
    else:
        ul_x, new_x_resolution, new_x_rot, ul_y, new_y_rot, new_y_resolution = hyObj.transform
        print(hyObj.projection,hyObj.map_info)

    transform_matrix = np.array([[new_x_resolution, new_x_rot],[new_y_rot, new_y_resolution]])

    if hyObj.map_info[0].startswith("Geographic"):

        if mapobj.zone is None: # Not defined, assume to be geographic from csv / point_df
            xy_coord_array = point_coord_df[['img_x','img_y']].values-np.array([[ul_x,ul_y]])
        else: # assume to has both UTM and coord in point_df
            xy_coord_array = point_coord_df[['lon','lat']].values-np.array([[ul_x,ul_y]])
    elif hyObj.map_info[0].startswith("UTM"):

        if point_epsg_code is None: # latlon in point, but utm in image
            img_zone = hyObj.map_info[7]+hyObj.map_info[8][0]
            img_mapobj = BasicMapObj(zone=img_zone) #NAD83_WGS84_obj,

            x_coord, y_coord = img_mapobj.convert_xycoord_gdal(point_coord_df['img_y'].values, point_coord_df['img_x'].values)
            xy_coord_array = np.stack((x_coord, y_coord)).T -np.array([[ul_x,ul_y]])
        else:
            xy_coord_array = point_coord_df[['img_x','img_y']].values-np.array([[ul_x,ul_y]])


    img_loc_array = (xy_coord_array@(np.linalg.inv(transform_matrix).T)).astype(np.int32)   # zero-based

    n_neighbor = max(0,n_neighbor)
    if n_neighbor>=0:

        if n_neighbor==0:
            offset_arr_col = np.array([[1,0]])
            offset_arr_row = np.array([[1,0]])    
            uid_list = np.repeat(point_coord_df[uid].values,1)
            new_uid_list = np.tile([f'_{x}' for x in range(1)],img_loc_array.shape[0])                    

        if n_neighbor== 4:

            offset_arr_col = np.array([[1,0],
                                        [1,0],
                                        [1,-1],
                                        [1,1],
                                        [1,0]])
            offset_arr_row = np.array([[1,0],
                                        [1,-1],
                                        [1,0],
                                        [1,0],
                                        [1,1]])

            uid_list = np.repeat(point_coord_df[uid].values,5)
            new_uid_list = np.tile([f'_{x}' for x in range(5)],img_loc_array.shape[0])

        if n_neighbor== 8:
            offset_arr_col = np.array([[1,0],
                                        [1,0],
                                        [1,-1],
                                        [1,1],
                                        [1,0],
                                        [1,-1],
                                        [1,1],
                                        [1,-1],
                                        [1,1]])
            offset_arr_row = np.array([[1,0],
                                        [1,-1],
                                        [1,0],
                                        [1,0],
                                        [1,1],
                                        [1,-1],
                                        [1,-1],
                                        [1,1],
                                        [1,1]])

            uid_list = np.repeat(point_coord_df[uid].values,9)
            new_uid_list = np.tile([f'_{x}' for x in range(9)],img_loc_array.shape[0])          

        img_loc_array_with_nb_col = offset_arr_col@np.vstack([img_loc_array[:,0],np.ones(img_loc_array.shape[0])])
        img_loc_array_with_nb_row = offset_arr_row@np.vstack([img_loc_array[:,1],np.ones(img_loc_array.shape[0])])
        new_uid_list = uid_list+new_uid_list

        img_loc_array_with_nb_col = img_loc_array_with_nb_col.T.ravel().astype(np.int32)
        img_loc_array_with_nb_row = img_loc_array_with_nb_row.T.ravel().astype(np.int32) # zero-based

        return_df = pd.DataFrame({'new_uid':new_uid_list,uid:uid_list,'img_col_glt':img_loc_array_with_nb_col,'img_row_glt':img_loc_array_with_nb_row})

    print('use_glt_bool',use_glt_bool)
    if use_glt_bool:
        valid_mask = (img_loc_array_with_nb_col>=0) & (img_loc_array_with_nb_col< hyObj.columns_glt) & (img_loc_array_with_nb_row>=0) & (img_loc_array_with_nb_row< hyObj.lines_glt)

        if valid_mask.sum()==0:
            print("No valid GLT locations.")
            return pd.DataFrame()

        return_df = return_df[valid_mask]

        post_glt_col_ind = hyObj.glt_x[(img_loc_array_with_nb_row[valid_mask],img_loc_array_with_nb_col[valid_mask])]-1
        post_glt_row_ind = hyObj.glt_y[(img_loc_array_with_nb_row[valid_mask],img_loc_array_with_nb_col[valid_mask])]-1 # one-based to zero-based

        return_df["img_col_raw"] = post_glt_col_ind.astype(np.int32)
        return_df["img_row_raw"] = post_glt_row_ind.astype(np.int32) # zero-based
    else:
        return_df["img_col_raw"] = return_df['img_col_glt']
        return_df["img_row_raw"] = return_df['img_row_glt']

    # check whether points are within the boundary of the image or not
    return_df = return_df[(return_df['img_col_raw']>=0) & (return_df['img_col_raw']< hyObj.columns) & (return_df['img_row_raw']>=0) & (return_df['img_row_raw']< hyObj.lines)]
    return return_df

def add_df_lat_lon(point_coord_neighbor_df, hyObj, mapobj, offset=0.5, use_glt_bool = False):
    ''' Add LAT LON of the points in the dataframe
  
    '''

    if use_glt_bool:
        ul_x, new_x_resolution, new_x_rot, ul_y, new_y_rot, new_y_resolution = hyObj.glt_transform
    else:
        ul_x, new_x_resolution, new_x_rot, ul_y, new_y_rot, new_y_resolution = hyObj.transform

    transform_matrix = np.array([[new_x_resolution, new_x_rot],[new_y_rot, new_y_resolution]])

    loc_array = point_coord_neighbor_df[['img_col_glt','img_row_glt']].values.transpose() # zero-based

    img_coord_array = np.dot(transform_matrix,loc_array+offset)+np.array([[ul_x],[ul_y]])

    if hyObj.map_info[0].startswith("Geographic"):
        point_coord_neighbor_df['lat'] = img_coord_array[1,:]
        point_coord_neighbor_df['lon'] = img_coord_array[0,:]
    elif hyObj.map_info[0].startswith("UTM"):
        lat_list,lon_list = mapobj.convert_latlon(img_coord_array[0,:],img_coord_array[1,:])
        point_coord_neighbor_df['lat'] = lat_list
        point_coord_neighbor_df['lon'] = lon_list


def subset_band_list(hyObj,spec_df,use_band_list, band_list):

    # do not subset bands, do nothing
    if use_band_list==False:
        return spec_df

    # subset bands
    else:
        # user does not provide band list, use bad band list as default
        if len(band_list)==0:
            # no bad band list in the file, do nothing
            if not isinstance(hyObj.bad_bands,np.ndarray):
                return spec_df
            # use bad band list
            else:
                return spec_df.iloc[:,hyObj.bad_bands]
        # user provides band list
        else:
            return spec_df.iloc[:, band_list]

def local_point2spec(hyObj, point_csv, uid, xcoord, ycoord, point_epsg_code, n_neighbor=4, use_band_list=True, band_list=[],use_glt_bool=False):
    """Extract spectra with points in a CSV from the hyperspectral image

    Parameters
    ----------
    hyObj :                   HyTools file object
    point_csv:              str
                                    full filename of the point CSV 
    uid:                        str
                                    the user specified unique point ID in the CSV
    xcoord:                  str
                                    the column name in CSV for X coordinate of the points
    ycoord:                  str
                                    the column name in CSV for Y coordinate of the points
    point_epsg_code:  int
                                    EPSG code for the projection of the points, XY coordinates are based on this projection
    n_neighbor:           int
                                    default is 4, other options are 0, 8
                                    how many neighbors in the image should be sampled from the center                           
    use_band_list:       boolean
                                    default True; whether to use a subset of bands                           
    band_list:               list or numpy array                         
                                    default is a blank list
                                    if it is a list, it should be one like [5,6,7,8,9, 12]
                                    if it is a numpy array, it should be the same size as hyObj.bad_bands with only True or False in the array

    use_glt_bool:      boolean
                                    default False; whether to use geo-lookup table for pixel indexing         
            
    Returns
    -------
    point_coord_neighbor_df: pandas dataframe
                                                it include all the location and spectra information for all points from the CSV
    
    """

    point_df = pd.read_csv(point_csv, sep=',')
    if point_epsg_code is None:
        if hyObj.map_info[0].startswith("UTM"):
            img_zone = hyObj.map_info[7]+hyObj.map_info[8][0]
            parameter_obj = BasicMapObj(zone=img_zone)  #NAD83_WGS84_obj,
        else:
            parameter_obj = BasicMapObj() #NAD83_WGS84_obj
    else:
        parameter_obj = BasicMapObj(zone=point_epsg_code) #NAD83_WGS84_obj,

    # create a dataframe with image georeferenced coordinates of all points of interest
    point_coord_df = local_transform_all_point(parameter_obj, point_df, uid, xcoord, ycoord,point_epsg_code)

    # create a dataframe with columns and lines of all image space neighbors of points of interest
    point_coord_neighbor_df = get_neighbor(hyObj, point_coord_df, n_neighbor, uid,point_epsg_code,parameter_obj,use_glt_bool)

    if point_coord_neighbor_df.shape[0]==0:
        print("0 point within boundary!\n\n")
        return None
    else:
        # add LAT LON of the points in the dataframe
        add_df_lat_lon(point_coord_neighbor_df, hyObj,parameter_obj,use_glt_bool=use_glt_bool)

        spec_data = hyObj.get_pixels(point_coord_neighbor_df['img_row_raw'].values,point_coord_neighbor_df['img_col_raw'].values) # zero-based

        # determine the column names of the spectra dataframe based on wavelengths
        if hyObj.wavelength_units.lower()[:4]=='micr':
            new_band_name = ['B{:0.3f}'.format(x) for x in hyObj.wavelengths]
        elif hyObj.wavelength_units.lower()[:4]=='nano' :
            new_band_name = ['B{:04d}'.format(int(x)) for x in hyObj.wavelengths]
        else:
            new_band_name = ['B{:d}'.format(x+1) for x in range(hyObj.bands)]

        spec_df = pd.DataFrame(spec_data, columns=new_band_name)

        # perform the subsetting of the columns in the dataframe according to the band_list or hyObj.bad_bands
        spec_df = subset_band_list(hyObj,spec_df,use_band_list, band_list)

        point_coord_neighbor_df=point_coord_neighbor_df.reset_index(drop=True)
        point_coord_neighbor_df = pd.concat([point_coord_neighbor_df,spec_df], axis=1, join='inner')

        return point_coord_neighbor_df

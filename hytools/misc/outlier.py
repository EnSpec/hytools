from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist  
import scipy.signal as sig




# Waves used for outlier/abnormal flight lines
waves = [480, 560, 660, 850, 950, 1050, 1240, 1650, 2217]

# Convert to band number
bands = [ray.get(actors[0].wave_to_band.remote(wave)) for wave in waves]

def bad_band_updater(hy_obj):
    hy_obj.bad_bands = ~hy_obj.bad_bands
    hy_obj.bad_bands[bands] = False

_ = ray.get([a.do.remote(bad_band_updater) for a in actors])

#For testing, delete
#bad_bands = ray.get(actors[0].do.remote(lambda x: x.bad_bands))

# Correct 
if correction =='topo':
    topo_coeffs(actors,config_dict['topo'])

group_dict = {}


for i in range(len(hyObj_pointer_dict_list)):
    print(hyObj_pointer_dict_list[i].file_name)
    if hyObj_pointer_dict_list[i].file_type == "ENVI":

        if hyObj_pointer_dict_list[i].interleave == 'bsq':
            spec_data = hyObj_pointer_dict_list[i].data[:, idxRand_dict[i][0], idxRand_dict[i][1]].transpose()
        elif hyObj_pointer_dict_list[i].interleave == 'bil':
            spec_data = hyObj_pointer_dict_list[i].data[idxRand_dict[i][0], :, idxRand_dict[i][1]]
        # hyObj.interleave=='bip':
        else:
            spec_data = hyObj_pointer_dict_list[i].data[idxRand_dict[i][0], idxRand_dict[i][1], :]
    elif hyObj_pointer_dict_list[i].file_type == "HDF":
        spec_data = hyObj_pointer_dict_list[i].data[idxRand_dict[i][0], idxRand_dict[i][1], :]
    else:
        return None

    wave_samples = spec_data[:, band_subset_outlier]
    wave_samples = wave_samples / image_smooth[i][band_subset_outlier]

    sub_index_img_tag = (sample_img_tag == i + 1)
    sample_cos_i_sub = sample_cos_i[sub_index_img_tag]
    sample_slope_sub = sample_slope[sub_index_img_tag]
    sample_c1_sub = sample_c1[sub_index_img_tag]

    topo_mask_sub = (sample_cos_i_sub > COSINE_I_MIN_THRESHOLD) & (sample_slope_sub > SLOPE_MIN_THRESHOLD)

    for iband in range(len(band_subset_outlier)):
        wave_samples_band = wave_samples[:, iband]

        topo_coeff, _, _ = generate_topo_coeff_band(wave_samples_band, (wave_samples_band > REFL_MIN_THRESHOLD) & (wave_samples_band < REFL_MAX_THRESHOLD) & topo_mask_sub, sample_cos_i_sub, non_negative=True)
        correctionFactor = (sample_c1_sub + topo_coeff) / (sample_cos_i_sub + topo_coeff)
        correctionFactor = correctionFactor * topo_mask_sub + 1.0 * (1 - topo_mask_sub)
        wave_samples[:, iband] = wave_samples_band * correctionFactor

    wave_all_samples = np.hstack((wave_all_samples, wave_samples.T))

ndvi_mask = (sample_ndvi > 0.15) & (sample_ndvi <= 0.95)
obs_mask = np.isfinite(sample_k_vol) & np.isfinite(sample_k_geom)
temp_mask = (wave_all_samples[0] > REFL_MIN_THRESHOLD) & (wave_all_samples[0] < REFL_MAX_THRESHOLD) & (obs_mask) & (ndvi_mask)

for iband in range(len(band_subset_outlier)):
    new_df = pd.DataFrame({'k_geom': sample_k_geom[temp_mask], 'k_vol': sample_k_vol[temp_mask],
                          'reflectance': wave_all_samples[iband, temp_mask], 'line_id': sample_img_tag[temp_mask],
                          "NDVI": sample_ndvi[temp_mask]})

    new_df['ndvi_cut_bins'] = pd.cut(new_df['NDVI'],
                                    bins=[0.15, 0.4, 0.7, 0.95],
                                    labels=['ndvi_1', 'ndvi_2', 'ndvi_3'])

    new_df['geom_cut_bins'] = pd.cut(new_df['k_geom'],
                            bins=np.percentile(sample_k_geom[temp_mask], [5, 33, 67, 95]),  # [5,33,67,95] #[5,25,50,75,95]
                            labels=['k_geom_1', 'k_geom_2', 'k_geom_3'])  # ,'k_geom_4'

    new_df['vol_cut_bins'] = pd.cut(new_df['k_vol'],
                            bins=np.percentile(sample_k_vol[temp_mask], [5, 33, 67, 95]),   # [5,25,50,75,95] # [5,33,67,95]
                            labels=['k_vol_1', 'k_vol_2', 'k_vol_3'])  # 'k_vol_4'

    new_df_bin_group_mean = new_df.groupby(['vol_cut_bins', 'geom_cut_bins', 'ndvi_cut_bins', 'line_id']).median()  # mean()

    new_df_bin_group_mean.reset_index(inplace=True)

    n_bin = new_df_bin_group_mean.shape[0] // len(hyObj_pointer_dict_list)

    ss = new_df_bin_group_mean["reflectance"].values

    bin_avg_array = np.reshape(ss, (n_bin, len(hyObj_pointer_dict_list)))

    bin_mean = np.nanmedian(bin_avg_array, axis=1)
    inds = np.where(np.isnan(bin_avg_array))

    # Place column means in the indices. Align the arrays using take
    bin_avg_array[inds] = np.take(bin_mean, inds[0])

    bin_avg_array = bin_avg_array / bin_mean[:, np.newaxis]

    bin_avg_array = bin_avg_array[~np.isnan(bin_avg_array[:, 0])]

    # Y = pdist(bin_avg_array.T, 'seuclidean', V=None)
    Y = pdist(bin_avg_array.T, 'euclidean', V=None)
    # Y = pdist(bin_avg_array.T, 'canberra')

    print(Y)

    return_dict = {}

    # H_s = hierarchy.single(Y)
    H_s = hierarchy.complete(Y)
    T_ = hierarchy.fcluster(H_s, 1.2, criterion='distance')
    print("Cluster thres 1.2", T_)

    return_dict["Cluster thres 1.2"] = T_.tolist()

    T_ = hierarchy.fcluster(H_s, 1.0, criterion='distance')
    print("Cluster thres 1.0", T_)

    return_dict["Cluster thres 1.0"] = T_.tolist()

    T_ = hierarchy.fcluster(H_s, 0.85, criterion='distance')
    print("Cluster thres 0.85", T_)

    return_dict["Cluster thres 0.9"] = T_.tolist()

    return_dict["distance of metrics"] = Y.tolist()

    major_label_id = np.bincount(np.array(T_)).argmax()

    outlier_img_tag = (np.array(T_) != major_label_id)

    return_dict["outlier_image_bool"] = outlier_img_tag.astype(int).tolist()
    return_dict["outlier_count"] = int(np.count_nonzero(outlier_img_tag))
    group_dict['b' + str(iband + 1)] = return_dict

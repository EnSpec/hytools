'''neon2envi.py

TODO: Add phase and UTC time to ancillary output
TODO: Implement progress bar like from this example:
      https://docs.ray.io/en/master/auto_examples/progress_bar.html

'''
import argparse
import pickle
import os
import ray
import numpy as np
from sklearn.decomposition import PCA
import hytools as ht
from hytools.io.envi import WriteENVI

def main():
    '''

    '''
    parser = argparse.ArgumentParser(description = "Perform a PCA")
    parser.add_argument('images',help="Input image pathnames", nargs='*')
    parser.add_argument('output_dir',help="Output directory", type = str)
    #parser.add_argument("-t", help="Transform type", type = str)
    parser.add_argument("-comps", help="Number of components to export", type = int,required=False,default=5)
    parser.add_argument("-sample", help="Percent of data to sample", type = float,required=False,default=0.1)
    parser.add_argument("-merge", help="Use gdal to mosaic", type = bool,required=False,default=False)

    args = parser.parse_args()

    if not args.output_dir.endswith("/"):
        args.output_dir+="/"

    if ray.is_initialized():
        ray.shutdown()
    ray.init(num_cpus = len(args.images))

    hytool = ray.remote(ht.HyTools)
    actors = [hytool.remote() for image in args.images]

    if args.images[0].endswith('.h5'):
        file_type = 'neon'
    else:
        file_type = 'envi'

    _ = ray.get([a.read_file.remote(image,file_type) for a,image in zip(actors,args.images)])

    # Sample data
    samples  = ray.get([a.do.remote(subsample,args) for a in actors])

    # Center, scale and fit PCA transform
    X = np.concatenate(samples).astype(np.float32)
    X /=X.mean(axis=1)[:,np.newaxis]
    X /=X.std(axis=1,ddof=1)[:,np.newaxis]
    X = X[~np.isnan(X.sum(axis=1)) & ~np.isinf(X.sum(axis=1)),:]

    print('Performing PCA decomposition')
    pca = PCA(n_components=args.comps)
    pca.fit(X)
    pca_pkl = pickle.dumps(pca)

    args.pca_pkl = pca_pkl
    #Apply tranform and export
    _  = ray.get([a.do.remote(apply_transform,args) for a in actors])

    if args.merge and len(args.images) > 1:
        print('Mosaicking flightlines')
        output_files = ["%s%s_pca" %(args.output_dir,image) for image in \
                        ray.get([a.do.remote(lambda x : x.base_name) for a in actors])]
        string = ['gdal_merge.py','-o', '%stranform_mosaic.tif' % output_dir] + output_files
        os.system(' '.join(string))

def subsample(hy_obj,args):

    print("Sampling %s" % os.path.basename(hy_obj.file_name))

    # Select 'sample_perc' % of pixels for modeling
    # This can probably be written more concisely
    sub_samples = np.zeros((hy_obj.lines,hy_obj.columns)).astype(bool)
    idx = np.array(np.where(hy_obj.mask['no_data'])).T
    idxRand= idx[np.random.choice(range(len(idx)),int(len(idx)*args.sample), replace = False)].T
    sub_samples[idxRand[0],idxRand[1]] = True
    hy_obj.mask['samples'] = sub_samples

    X = []

    hy_obj.create_bad_bands([[300,400],[1337,1430],[1800,1960],[2450,2600]])
    for band_num,band in enumerate(hy_obj.bad_bands):
        if ~band:
            X.append(hy_obj.get_band(band_num,mask='samples'))
    return  np.array(X).T

def apply_transform(hy_obj,args):

    print("Exporting %s PCA" % hy_obj.base_name)
    pca = pickle.loads(args.pca_pkl)
    output_name = '%s/%s_pca' % (args.output_dir,hy_obj.base_name)
    header_dict = hy_obj.get_header()
    header_dict['bands'] = pca.n_components
    header_dict['wavelength'] = []
    header_dict['fwhm'] = []
    header_dict['data type'] = 4
    header_dict['data ignore value'] = 0
    writer = WriteENVI(output_name,header_dict)
    iterator = hy_obj.iterate(by = 'chunk',chunk_size = (500,500))

    while not iterator.complete:
        chunk = iterator.read_next()

        X_chunk = chunk[:,:,~hy_obj.bad_bands].astype(np.float32)
        X_chunk = X_chunk.reshape((X_chunk.shape[0]*X_chunk.shape[1],X_chunk.shape[2]))
        X_chunk /=X_chunk.mean(axis=1)[:,np.newaxis]
        X_chunk /= X_chunk.std(axis=1,ddof=1)[:,np.newaxis]
        X_chunk[np.isnan(X_chunk) | np.isinf(X_chunk)] = 0

        pca_chunk=  pca.transform(X_chunk)
        pca_chunk = pca_chunk.reshape((chunk.shape[0],chunk.shape[1],pca.n_components))
        pca_chunk[chunk[:,:,0] == hy_obj.no_data] =0

        writer.write_chunk(pca_chunk,
                           iterator.current_line,
                           iterator.current_column)






if __name__== "__main__":
    main()

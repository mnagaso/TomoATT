# functinos for post-processing output data in h5 file

import numpy
import h5py


# output arrays has one overlapping layer between adjusted subdomains
# this function will skip those layers and reconstruct the entire grid dataset
def get_data_from_h5(fpath, dataset_name, nr_glob, nt_glob, np_glob, ndiv_r, ndiv_t, ndiv_p, verbose=False):
    """

    fpath: path to h5 file
    dataset_name: name of dataset in h5 file
    nr_glob: number of grid points in r direction
    nt_glob: number of grid points in t direction
    np_glob: number of grid points in p direction
    ndiv_r: number of subdomains in r direction
    ndiv_t: number of subdomains in t direction
    ndiv_p: number of subdomains in p direction
    verbose: print out information

    """
    # total number of subdomains
    n_sub = ndiv_r*ndiv_t*ndiv_p

    # total points on each direction with overlap
    nr_total_glob = nr_glob + ndiv_r - 1
    nt_total_glob = nt_glob + ndiv_t - 1
    np_total_glob = np_glob + ndiv_p - 1

    # prepare a 3D array to store the data
    data_glob = numpy.zeros((nr_glob,nt_glob,np_glob), dtype=numpy.float64)

    # read h5 file
    with h5py.File(fpath, 'r') as f:
        # load data data by each subdomain

        # offset
        offset = 0

        for ir_sub in range(ndiv_r):
            for it_sub in range(ndiv_t):
                for ip_sub in range(ndiv_p):

                    # number of data point for this sub domain
                    nr_sub = nr_glob//ndiv_r
                    nt_sub = nt_glob//ndiv_t
                    np_sub = np_glob//ndiv_p

                    # offset for each direction
                    offset_r = ir_sub*nr_sub
                    offset_t = it_sub*nt_sub
                    offset_p = ip_sub*np_sub

                    # add modulus to the last subdomains
                    if ir_sub == ndiv_r-1:
                        nr_sub += nr_glob%ndiv_r
                    if it_sub == ndiv_t-1:
                        nt_sub += nt_glob%ndiv_t
                    if ip_sub == ndiv_p-1:
                        np_sub += np_glob%ndiv_p

                    # add overlap layer if this subdomain is not the last one for each direction
                    if ir_sub != ndiv_r-1:
                        nr_sub += 1
                    if it_sub != ndiv_t-1:
                        nt_sub += 1
                    if ip_sub != ndiv_p-1:
                        np_sub += 1

                    # number of data point for this sub domain
                    n_points_total_sub = nr_sub*nt_sub*np_sub

                    # load data
                    data_sub = f[dataset_name][offset:offset+n_points_total_sub]

                    # reshape data
                    data_sub = data_sub.reshape(nr_sub, nt_sub, np_sub)

                    # put those data in global 3d array
                    data_glob[offset_r:offset_r+nr_sub, offset_t:offset_t+nt_sub, offset_p:offset_p+np_sub] = data_sub

                    # update offset
                    offset += n_points_total_sub


    return data_glob


# output arrays has one overlapping layer between adjusted subdomains
# this function will skip those layers and reconstruct the entire grid dataset
def get_data_from_ascii(fpath, nr_glob, nt_glob, np_glob, ndiv_r, ndiv_t, ndiv_p, verbose=False):
    """

    fpath: path to ascii data file
    nr_glob: number of grid points in r direction
    nt_glob: number of grid points in t direction
    np_glob: number of grid points in p direction
    ndiv_r: number of subdomains in r direction
    ndiv_t: number of subdomains in t direction
    ndiv_p: number of subdomains in p direction
    verbose: print out information

    """
    # total number of subdomains
    n_sub = ndiv_r*ndiv_t*ndiv_p

    # total points on each direction with overlap
    nr_total_glob = nr_glob + ndiv_r - 1
    nt_total_glob = nt_glob + ndiv_t - 1
    np_total_glob = np_glob + ndiv_p - 1

    # prepare a 3D array to store the data
    data_glob = numpy.zeros((nr_glob,nt_glob,np_glob), dtype=numpy.float64)

    # read data
    data_tmp = numpy.loadtxt(fpath)

    # load data data by each subdomain

    # offset
    offset = 0

    for ir_sub in range(ndiv_r):
        for it_sub in range(ndiv_t):
            for ip_sub in range(ndiv_p):

                # number of data point for this sub domain
                nr_sub = nr_glob//ndiv_r
                nt_sub = nt_glob//ndiv_t
                np_sub = np_glob//ndiv_p

                # offset for each direction
                offset_r = ir_sub*nr_sub
                offset_t = it_sub*nt_sub
                offset_p = ip_sub*np_sub

                # add modulus to the last subdomains
                if ir_sub == ndiv_r-1:
                    nr_sub += nr_glob%ndiv_r
                if it_sub == ndiv_t-1:
                    nt_sub += nt_glob%ndiv_t
                if ip_sub == ndiv_p-1:
                    np_sub += np_glob%ndiv_p

                # add overlap layer if this subdomain is not the last one for each direction
                if ir_sub != ndiv_r-1:
                    nr_sub += 1
                if it_sub != ndiv_t-1:
                    nt_sub += 1
                if ip_sub != ndiv_p-1:
                    np_sub += 1

                # number of data point for this sub domain
                n_points_total_sub = nr_sub*nt_sub*np_sub

                # load data
                data_sub = data_tmp[offset:offset+n_points_total_sub]

                # reshape data
                data_sub = data_sub.reshape(nr_sub, nt_sub, np_sub)

                # put those data in global 3d array
                data_glob[offset_r:offset_r+nr_sub, offset_t:offset_t+nt_sub, offset_p:offset_p+np_sub] = data_sub

                # update offset
                offset += n_points_total_sub


    return data_glob






if __name__ == '__main__':

    # examples
    fpath = '../output/tomoatt_output_data_test.h5'
    nr = 10
    nt = 10
    np = 10
    ndiv_r = 2
    ndiv_t = 2
    ndiv_p = 2
    Ks_tomoatt = get_data_from_h5(fpath, "Data/Ks_inv_0000", nr, nt, np, ndiv_r, ndiv_t, ndiv_p, verbose=True)
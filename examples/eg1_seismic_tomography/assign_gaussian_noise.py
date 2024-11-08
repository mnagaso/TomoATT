# %%
# from pytomoatt.src_rec import SrcRec

# def assign_noise_to_src_rec_file(in_fname, out_fname, noise_level=0.1):
#     sr = SrcRec.read(in_fname)
#     sr.add_noise(noise_level)
#     sr.write(out_fname)


# if __name__ == "__main__":
#     in_fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward.dat" # input source receiver file
#     out_fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward_noisy.dat" # output source receiver file
#     sigma = 0.1 # noise level in seconds
#     assign_noise_to_src_rec_file(in_fname, out_fname, noise_level=sigma)

# %%
import sys
sys.path.append('../utils')
import functions_for_data as ffd

# read src_rec_file
ev, st = ffd.read_src_rec_file('OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward.dat')

# assign gaussian noise to the data
# ffd.assign_gaussian_noise(ev, 0.1)

# output file
ffd.write_src_rec_file('OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward_noisy.dat', ev, st)

# # perturbate earthquake locations
# lat_pert = 0.1 # degrees
# lon_pert = 0.1 # degrees
# dep_pert = 10  # km
# tau_pert = 0.5 # seconds
# ffd.assign_uniform_noise_to_ev(ev, lat_pert, lon_pert, dep_pert, tau_pert)

# # output file
# ffd.write_src_rec_file('OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward_errloc.dat', ev, st)



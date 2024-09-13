import sys
sys.path.append('../utils')
import functions_for_data as tools

fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward.dat"
ev_info, st_info = tools.read_src_rec_file(fname)

sigma = 0.1
ev_info = tools.assign_gaussian_noise(ev_info,sigma)

fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward_noisy.dat"
tools.write_src_rec_file(fname, ev_info, st_info)
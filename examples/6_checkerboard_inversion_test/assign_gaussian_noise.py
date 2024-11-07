from pytomoatt.src_rec import SrcRec


def assign_noise_to_src_rec_file(in_fname, out_fname, noise_level=0.1):
    sr = SrcRec.read(in_fname)
    sr.add_noise(noise_level)
    sr.write(out_fname)


if __name__ == "__main__":
    in_fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward.dat" # input source receiver file
    out_fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward_noisy.dat" # output source receiver file
    sigma = 0.1 # noise level in seconds
    assign_noise_to_src_rec_file(in_fname, out_fname, noise_level=sigma)
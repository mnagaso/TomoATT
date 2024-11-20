from pytomoatt.src_rec import SrcRec

class AssignNoise:
    def __init__(self, in_fname, out_fname):
        self.in_fname = in_fname
        self.out_fname = out_fname
        self.sr = SrcRec.read(self.in_fname)

    def assign_noise_for_tt(self, noise_level=0.1):
        self.sr.add_noise(noise_level)

    def assign_noise_for_src(self, lat_pert=0.1, lon_pert=0.1, dep_pert=10, tau_pert=0.5):
        self.sr.add_noise_to_source(lat_pert, lon_pert, dep_pert, tau_pert)

if __name__ == "__main__":
    in_fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward.dat" # input source receiver file
    out_fname = "OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward_errloc.dat" # output source receiver file
    sigma = 0.1 # noise level in seconds
    lat_pert = 0.1 # assign noise for latitude in degrees
    lon_pert = 0.1 # assign noise for longitude in degrees
    dep_pert = 10 # assign noise for depth in km
    tau_pert = 0.5 # assign noise for origin time in seconds

    # Initialize the instance
    an = AssignNoise(in_fname, out_fname)

    # Assign noise for travel time
    an.assign_noise_for_tt(sigma)
    
    # Assign noise for source
    an.assign_noise_for_src(lat_pert, lon_pert, dep_pert, tau_pert)
    
    # Write the output file
    an.sr.write(out_fname)
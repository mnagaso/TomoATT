import numpy as np
import os
import h5py
from pytomoatt.para import ATTPara


class BuildInitialModel():
    def __init__(self, par_file="./3_input_params/0_input_params_forward_simulation.yaml", output_dir="2_models"):
        """
        Build initial model for tomography inversion
        """
        self.ip = ATTPara(par_file)
        self.dep, self.lat, self.lon, self.dd, self.dt, self.dp = self.ip.init_axis()
        self.output_dir = output_dir

    def build_initial_model(self):
        """
        Build initial model for tomography inversion
        """
        n_rtp = self.ip.input_params['domain']['n_rtp']
        self.vel_init    = np.zeros(n_rtp)   # velocity
        self.xi_init     = np.zeros(n_rtp)   # xi and eta are azimuthal anisotropy
        self.eta_init    = np.zeros(n_rtp)
        self.zeta_init   = np.zeros(n_rtp)   # zeta is radial anisotropy

        for ir in range(n_rtp[0]):
            for it in range(n_rtp[1]):
                for ip in range(n_rtp[2]):
                    dep = self.dep[ir]
                    if (dep < 0):
                        self.vel_init[ir,it,ip]  = 6.0
                    elif (dep >= 0  and dep < 40):
                        self.vel_init[ir,it,ip]  = 6.0 + dep/40*2.0
                    else:
                        self.vel_init[ir,it,ip]  = 8.0

                    self.xi_init[ir,it,ip]   = 0.0
                    self.eta_init[ir,it,ip]  = 0.0
                    self.zeta_init[ir,it,ip] = 0.0


    def write(self):
        """
        write out in hdf5 format
        """
        os.makedirs(self.output_dir, exist_ok=True)
        with h5py.File('{}/model_init_N{:d}_{:d}_{:d}.h5'.format(self.output_dir, *self.ip.input_params['domain']['n_rtp']), 'w') as fout_init:
        # write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
            fout_init.create_dataset('eta', data=self.eta_init)
            fout_init.create_dataset('xi', data=self.xi_init)
            fout_init.create_dataset('zeta', data=self.zeta_init)
            fout_init.create_dataset('vel', data=self.vel_init)



if __name__ == "__main__":
    bim = BuildInitialModel()
    bim.build_initial_model()
    bim.write()

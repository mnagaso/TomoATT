import numpy as np
import os
import h5py
from pytomoatt.para import ATTPara


class BuildCkbModel():
    def __init__(self, par_file="./3_input_params/0_input_params_forward_simulation.yaml", output_dir="2_models"):
        """
        Build ckb model for tomography inversion
        """
        self.ip = ATTPara(par_file)
        self.dep, self.lat, self.lon, self.dd, self.dt, self.dp = self.ip.init_axis()
        self.output_dir = output_dir

    def build_ckb_model(self):
        """
        Build ckb model for tomography inversion
        """
        vel_pert = 0.05
        ani_pert = 0.05

        # read initial model
        model_fname = "2_models/model_init_N61_61_61.h5"   # 检测板反演的最终模型
        with h5py.File(model_fname, 'r') as fmodel:
            vel_init    = np.array(fmodel['vel'])


        n_rtp = self.ip.input_params['domain']['n_rtp']
        self.vel_ckb    = np.zeros(n_rtp)   # velocity
        self.xi_ckb     = np.zeros(n_rtp)   # xi and eta are azimuthal anisotropy
        self.eta_ckb    = np.zeros(n_rtp)
        self.zeta_ckb   = np.zeros(n_rtp)   # zeta is radial anisotropy

        for ir in range(n_rtp[0]):
            for it in range(n_rtp[1]):
                for ip in range(n_rtp[2]):
                    dep = self.dep[ir]
                    lat = self.lat[it]
                    lon = self.lon[ip]

                    if (dep >= 0  and dep <= 40):
                        sigma   = np.sin(np.pi*lat/(0.5)) * np.sin(np.pi*lon/(0.5)) * np.sin(np.pi*dep/(20))
                        if sigma < 0:
                            psi = 60.0/180.0*np.pi
                        elif sigma > 0:
                            psi = 150.0/180.0*np.pi
                    else:
                        sigma = 0.0
                        psi   = 0.0

                    self.vel_ckb[ir,it,ip]   = vel_init[ir,it,ip] * (1.0 + vel_pert * sigma)
                    self.xi_ckb[ir,it,ip]    = ani_pert * abs(sigma) * np.cos(2*psi) 
                    self.eta_ckb[ir,it,ip]   = ani_pert * abs(sigma) * np.sin(2*psi) 
                    self.zeta_ckb[ir,it,ip]  = 0.0



    def write(self):
        """
        write out in hdf5 format
        """
        os.makedirs(self.output_dir, exist_ok=True)
        with h5py.File('{}/model_ckb_N{:d}_{:d}_{:d}.h5'.format(self.output_dir, *self.ip.input_params['domain']['n_rtp']), 'w') as fout_init:
        # write out the arrays eta_init, xi_init, zeta_init, fun_init, a_init, b_init, c_init, f_init
            fout_init.create_dataset('eta', data=self.eta_ckb)
            fout_init.create_dataset('xi', data=self.xi_ckb)
            fout_init.create_dataset('zeta', data=self.zeta_ckb)
            fout_init.create_dataset('vel', data=self.vel_ckb)



if __name__ == "__main__":
    bim = BuildCkbModel()
    bim.build_ckb_model()
    bim.write()

# %%
# download src_ref_files from Zenodo
import os
import requests

url = 'https://zenodo.org/records/14053821/files/src_rec_config.dat?download=1'

path = "1_src_rec_files/src_rec_config.dat"

# check file existence
if not os.path.exists(path):
    try:
        os.mkdir("1_src_rec_files")
    except:
        pass
    print("Downloading src_rec_config.dat from Zenodo...")
    response = requests.get(url, stream=True)
    with open(path, 'wb') as out_file:
        out_file.write(response.content)
    print("Download complete.")
else:
    print("src_rec_config.dat already exists.")

# %%
# prepare initial model

import numpy as np
import os
import h5py
from pytomoatt.model import ATTModel


class BuildInitialModel():
    def __init__(self, par_file="./3_input_params/input_params_signal.yaml", output_dir="2_models"):
        """
        Build initial model for tomography inversion
        """
        self.am = ATTModel(par_file)
        self.output_dir = output_dir

    def build_initial_model(self):
        """
        Build initial model for tomography inversion
        """
        for ir, dep in enumerate(self.am.depths):
            if (dep < 0):
                self.am.vel[ir, :, :] = 5.0
            elif (0 <= dep < 40.0):
                self.am.vel[ir, :, :] = 5.0 + dep/40*3.0
            else:
                self.am.vel[ir, :, :] = 8.0

if __name__ == "__main__":
    try:
        os.mkdir("2_models")
    except:
        pass
    bim = BuildInitialModel()
    bim.build_initial_model()
    bim.am.write('{}/model_init_N{:d}_{:d}_{:d}.h5'.format(bim.output_dir, *bim.am.n_rtp))

# %%
class BuildCkbModel():
    def __init__(self, par_file="./3_input_params/input_params_signal.yaml", output_dir="2_models"):
        """
        Build initial model for tomography inversion
        """
        self.am = ATTModel(par_file)
        self.output_dir = output_dir

    def build_ckb_model(self):
        """
        Build ckb model for tomography inversion
        """
        vel_pert = 0.2
        ani_pert = 0.1

        # read initial model
        model_fname = "2_models/model_init_N61_61_61.h5"   # 检测板反演的最终模型
        with h5py.File(model_fname, 'r') as fmodel:
            vel_init    = np.array(fmodel['vel'])

        self.am.vel = vel_init

        for ir, dep in enumerate(self.am.depths):
            for ip, lat in enumerate(self.am.latitudes):
                for it, lon in enumerate(self.am.longitudes):
                    if (dep >= 0  and dep <= 40 and lat >= 0.5 and lat <= 1.5 and lon >= 0.5 and lon <= 1.5):
                        sigma   = np.sin(np.pi*lat/(0.5)) * np.sin(np.pi*lon/(0.5)) * np.sin(np.pi*dep/(40))
                        if sigma < 0:
                            psi = 60.0/180.0*np.pi
                        elif sigma > 0:
                            psi = 150.0/180.0*np.pi
                    else:
                        sigma = 0.0
                        psi   = 0.0

                    self.am.vel[ir,it,ip]   = vel_init[ir,it,ip] * (1.0 + vel_pert * sigma)
                    self.am.xi[ir,it,ip]    = ani_pert * abs(sigma) * np.cos(2*psi) 
                    self.am.eta[ir,it,ip]   = ani_pert * abs(sigma) * np.sin(2*psi) 
                    self.am.zeta[ir,it,ip]  = 0.0


if __name__ == "__main__":
    try:
        os.mkdir("2_models")
    except:
        pass
    bim = BuildCkbModel()
    bim.build_ckb_model()
    bim.am.write('{}/model_ckb_N{:d}_{:d}_{:d}.h5'.format(bim.output_dir, *bim.am.n_rtp))




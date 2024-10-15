import numpy as np
import os
import h5py
from pytomoatt.model import ATTModel


class BuildInitialModel():
    def __init__(self, par_file="./3_input_params/0_input_params_forward_simulation.yaml", output_dir="2_models"):
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
                self.am.vel[ir, :, :] = 6.0
            elif (0 <= dep < 40.0):
                self.am.vel[ir, :, :] = 6.0 + dep/40*2.0
            else:
                self.am.vel[ir, :, :] = 8.0


if __name__ == "__main__":
    os.mkdir("2_models")
    bim = BuildInitialModel()
    bim.build_initial_model()
    bim.am.write('{}/model_init_N{:d}_{:d}_{:d}.h5'.format(bim.output_dir, *bim.am.n_rtp))

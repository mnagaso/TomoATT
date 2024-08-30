import numpy as np
import os
import h5py
from pytomoatt.para import ATTPara
from pytomoatt.checkerboard import Checker


class BuildCheckerboardModel():
    def __init__(self, par_file="./3_input_params/0_input_params_forward_simulation.yaml", output_dir="2_models") -> None:
        self.ip = ATTPara(par_file)
        self.dep, self.lat, self.lon, self.dd, self.dt, self.dp = self.ip.init_axis()
        self.output_dir = output_dir
    
    def build_checkerboard_model(self, npert=[2,2,2], vel_pert=0.06, ani_pert=0.06, margin=0.5, lim_z=[0, 40]):
        """
        Build checkerboard model for tomography inversion
        
        :param npert: Number of perturbations in each direction, defaults to [2,2,2]
        :type npert: list, optional
        :param vel_pert: Velocity perturbation, defaults to 0.06
        :type vel_pert: float, optional
        :param ani_pert: Anisotropy perturbation, defaults to 0.06
        :type ani_pert: float, optional
        :param margin: Margin along longitude and latitude direction, defaults to 0.5
        :type margin: float, optional
        :param lim_z: Depth limit, defaults to [0, 40]
        :type lim_z: list, optional
        """
        self.checker = Checker('./2_models/model_init_N61_61_61.h5')
        self.checker.init_axis(
            self.ip.input_params['domain']['min_max_dep'],
            self.ip.input_params['domain']['min_max_lat'],
            self.ip.input_params['domain']['min_max_lon'],
            self.ip.input_params['domain']['n_rtp'],
        )
        self.checker.checkerboard(
            *npert,
            pert_vel=vel_pert, pert_ani=ani_pert,
            lim_x=[self.ip.input_params['domain']['min_max_lon'][0]+margin, self.ip.input_params['domain']['min_max_lon'][1]-margin],
            lim_y=[self.ip.input_params['domain']['min_max_lat'][0]+margin, self.ip.input_params['domain']['min_max_lat'][1]-margin],
            lim_z=lim_z
        )
        self.checker.write("{}/model_ckb_N{:d}_{:d}_{:d}.h5".format(self.output_dir, *self.ip.input_params['domain']['n_rtp']))

if __name__ == "__main__":
    bcm = BuildCheckerboardModel()
    bcm.build_checkerboard_model()
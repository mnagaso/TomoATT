import numpy as np
from datetime import datetime
import pandas as pd
from pytomoatt.src_rec import SrcRec
from pytomoatt.para import ATTPara
import os

OTIME = datetime(1970, 1, 1, 0, 0, 0)

class GenerateSrcRecFile():
    def __init__(self, nsrc_x, nsrc_y:int, nrec_x:int, nrec_y:int,
                 par_file="./3_input_params/0_input_params_forward_simulation.yaml",
                 output_dir="./1_src_rec_files") -> None:
        self.ip = ATTPara(par_file)
        self.sr = SrcRec('')
        self.output_dir = output_dir
        self.nsrc_x = nsrc_x
        self.nsrc_y = nsrc_y
        self.nrec_x = nrec_x
        self.nrec_y = nrec_y

    def generate_src_rec(self, dep_range=[2, 35]):
        srcs = []
        count = 0
        for i in range(self.nsrc_x):
            for j in range(self.nsrc_y):
                for k in range(3):
                    src = {
                        "origin_time": OTIME,
                        "evla": self.ip.input_params['domain']['min_max_lat'][0] + 0.2 + i/(self.nsrc_x-1)*1.6,
                        "evlo": self.ip.input_params['domain']['min_max_lat'][0] + 0.2 + j/(self.nsrc_y-1)*1.6,
                        "evdp": 10 + k*10,
                        "mag": 2.0,
                        "num_rec": 0,
                        "event_id": f"MX{count:04d}",
                        "weight": 1.0,
                    }
                    srcs.append(src)
                    count += 1
        self.sr.src_points = pd.DataFrame(srcs)

        recs = []
        for i in range(self.nrec_x):
            for j in range(self.nrec_y):
                rec = {
                    "staname": f"JC{i*self.nrec_y+j:02d}",
                    "stla": self.ip.input_params['domain']['min_max_lat'][0] + 0.2 + i/(self.nrec_x-1)*1.6,
                    "stlo": self.ip.input_params['domain']['min_max_lon'][0] + 0.2 + j/(self.nrec_y-1)*1.6,
                    "stel": 0.0,
                }
                recs.append(rec)
        self.sr.receivers = pd.DataFrame(recs)
        self.sr.receivers.reset_index(drop=True, inplace=True)

        # generate source receiver pairs
        data = []
        for i, src in self.sr.src_points.iterrows():
            count = 0
            for j, rec in self.sr.receivers.iterrows():
                trace = {
                    "src_index": i,
                    "rec_index": j,
                    "staname": rec['staname'],
                    "stla": rec['stla'],
                    "stlo": rec['stlo'],
                    "stel": rec['stel'],
                    "phase": "P",
                    "tt": 1.0,
                    "weight": 1.0,
                }
                count += 1
                data.append(trace)
            src['num_rec'] = count
        self.sr.rec_points = pd.DataFrame(data)

    def generate_double_diff_data(self, azimuth_gap=15, dist_gap_cs=0.5, dist_gap_cr=0.1):
        self.sr.generate_double_difference('cs', max_azi_gap=azimuth_gap, max_dist_gap=dist_gap_cs)
        self.sr.generate_double_difference('cr', max_azi_gap=azimuth_gap, max_dist_gap=dist_gap_cr)

    def write(self):
        os.makedirs(self.output_dir, exist_ok=True)
        self.sr.write(f"{self.output_dir}/src_rec_config.dat")


if __name__ == "__main__":
    gsr = GenerateSrcRecFile(nsrc_x=17, nsrc_y=17, nrec_x=5, nrec_y=5)
    gsr.generate_src_rec()
    gsr.generate_double_diff_data()
    gsr.write()
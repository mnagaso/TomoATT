# %%
import math
import numpy as np

RAD2DEG = 180/math.pi
DEG2RAD = math.pi/180
R_earth = 6371.0

# read src_rec_file (only for abs time)
def read_src_rec_file(filename):
    doc = open(filename,'r')
    doc_input = doc.readlines()
    doc.close()

    time = []
    ev_time = []
    cc = 0
    for info in doc_input:
        tmp=info.split()
        if(cc == 0):    # ev info line            
            ndata = int(tmp[11])
            cc = cc + 1
            tmp_ev_time = []
        else:
            if (len(tmp) < 10):     # abs data
                time.append(float(tmp[7]))
                tmp_ev_time.append(float(tmp[7]))
                if(cc == ndata):
                    ev_time.append(tmp_ev_time)
                    cc = 0
                else:
                    cc = cc + 1
            else:       # cs_dif or cr_dif
                time.append(float(tmp[12]))
                tmp_ev_time.append(float(tmp[12]))
                if(cc == ndata):
                    ev_time.append(tmp_ev_time)
                    cc = 0
                else:
                    cc = cc + 1
    return np.array(time),ev_time


# %%
fname = 'src_rec_obs.dat'
time_true,ev_time_true = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_no_swap_serial/src_rec_file_step_0000.dat'
time_no_swap_serial,ev_time_no_swap_serial = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_no_swap_parallel/src_rec_file_step_0000.dat'
time_no_swap_parallel,ev_time_no_swap_parallel = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_swap_serial/src_rec_file_step_0000.dat'
time_swap_serial,ev_time_swap_serial = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_swap_parallel/src_rec_file_step_0000.dat'
time_swap_parallel,ev_time_swap_parallel = read_src_rec_file(fname)



# %%
error_no_swap_serial = ((time_no_swap_serial-time_true))
if (error_no_swap_serial.std() < 0.1):
    print("For no_swap_serial,   src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))
else:
    print("For no_swap_serial,   output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_no_swap_parallel = ((time_no_swap_parallel-time_true))
if (error_no_swap_parallel.std() < 0.1):
    print("For no_swap_parallel, src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_parallel.std())))
else:
    print("For no_swap_parallel, output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))


error_swap_serial = ((time_swap_serial-time_true))
if (error_swap_serial.std() < 0.1):
    print("For swap_serial,      src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_serial.std())))
else:
    print("For swap_serial,      output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_swap_parallel = ((time_swap_parallel-time_true))
if (error_swap_parallel.std() < 0.1):
    print("For swap_parallel,    src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_parallel.std())))
else:
    print("For swap_parallel,    output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))




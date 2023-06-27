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
            time.append(float(tmp[7]))
            tmp_ev_time.append(float(tmp[7]))
            if(cc == ndata):
                ev_time.append(tmp_ev_time)
                cc = 0
            else:
                cc = cc + 1
    return np.array(time),ev_time


# %%
fname = 'src_rec_true.dat'
time_true,ev_time_true = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_N41_41_41/src_rec_file_forward.dat'
time_N41,ev_time_N41 = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_N61_61_61/src_rec_file_forward.dat'
time_N61,ev_time_N61 = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_N81_81_81/src_rec_file_forward.dat'
time_N81,ev_time_N81 = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_N121_121_121/src_rec_file_forward.dat'
time_N121,ev_time_N121 = read_src_rec_file(fname)

fname = 'OUTPUT_FILES_N161_161_161/src_rec_file_forward.dat'
time_N161,ev_time_N161 = read_src_rec_file(fname)


# %%
error_N41 = abs((time_N41-time_true))
print("The numerical error in a mesh of 41^3 is: %6.3f s (mean), %6.3f s (std), %6.3f (max)"%(error_N41.mean(),error_N41.std(),error_N41.max()))

error_N61 = abs((time_N61-time_true))
order = math.log(error_N61.mean()/error_N41.mean())/math.log(40/60)
print("The numerical error in a mesh of 61^3 is: %6.3f s (mean), %6.3f s (std), %6.3f (max), order of accuracy is %5.2f"%
      (error_N61.mean(),error_N61.std(),error_N61.max(),order))

error_N81 = abs((time_N81-time_true))
order = math.log(error_N81.mean()/error_N61.mean())/math.log(60/80)
print("The numerical error in a mesh of 81^3 is: %6.3f s (mean), %6.3f s (std), %6.3f (max), order of accuracy is %5.2f"%
      (error_N81.mean(),error_N81.std(),error_N81.max(),order))

error_N121 = abs((time_N121-time_true))
order = math.log(error_N121.mean()/error_N81.mean())/math.log(80/120)
print("The numerical error in a mesh of 121^3 is: %6.3f s (mean), %6.3f s (std), %6.3f (max), order of accuracy is %5.2f"%
      (error_N121.mean(),error_N121.std(),error_N121.max(),order))

error_N161 = abs((time_N161-time_true))
order = math.log(error_N161.mean()/error_N121.mean())/math.log(120/160)
print("The numerical error in a mesh of 161^3 is: %6.3f s (mean), %6.3f s (std), %6.3f (max), order of accuracy is %5.2f"%
      (error_N161.mean(),error_N161.std(),error_N161.max(),order))





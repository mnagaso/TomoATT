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


def check_result_times(errors, ngrid, order=None, tol=0.02):
    print("The numerical error in a mesh of %d^3 is: %6.3f s (mean), %6.3f s (std), %6.3f (max)" %
         (ngrid, errors.mean(), errors.std(), errors.max()))

    if order is not None:
        print("The order of accuracy is: %6.3f" % order)

    # check the order of accuracy is within tolerance
    if abs(order-1.0) < tol:
        print("The order of accuracy is within tolerance")
    else:
        print("The order of accuracy is NOT within tolerance")
        exit(1) # return non-zero exit code (test failed)

if __name__ == '__main__':

    # read theoretical travel times
    fname = 'src_rec_true.dat'
    time_true, ev_time_true = read_src_rec_file(fname)

    # read calculated travel times
    fname = 'OUTPUT_FILES_N41_41_41/src_rec_file_forward.dat'
    time_N41, ev_time_N41 = read_src_rec_file(fname)

    fname = 'OUTPUT_FILES_N61_61_61/src_rec_file_forward.dat'
    time_N61, ev_time_N61 = read_src_rec_file(fname)

    fname = 'OUTPUT_FILES_N81_81_81/src_rec_file_forward.dat'
    time_N81, ev_time_N81 = read_src_rec_file(fname)

    fname = 'OUTPUT_FILES_N121_121_121/src_rec_file_forward.dat'
    time_N121, ev_time_N121 = read_src_rec_file(fname)

    fname = 'OUTPUT_FILES_N161_161_161/src_rec_file_forward.dat'
    time_N161, ev_time_N161 = read_src_rec_file(fname)

    # calculate errors
    error_N41 = abs((time_N41-time_true))
    check_result_times(error_N41, 41)

    error_N61 = abs((time_N61-time_true))
    order = math.log(error_N61.mean()/error_N41.mean())/math.log(40/60)
    check_result_times(error_N61, 61, order)

    error_N81 = abs((time_N81-time_true))
    order = math.log(error_N81.mean()/error_N61.mean())/math.log(60/80)
    check_result_times(error_N81, 81, order)

    error_N121 = abs((time_N121-time_true))
    order = math.log(error_N121.mean()/error_N81.mean())/math.log(80/120)
    check_result_times(error_N121, 121, order)

    error_N161 = abs((time_N161-time_true))
    order = math.log(error_N161.mean()/error_N121.mean())/math.log(120/160)
    check_result_times(error_N161, 161, order)

    exit(0) # return zero exit code (success)




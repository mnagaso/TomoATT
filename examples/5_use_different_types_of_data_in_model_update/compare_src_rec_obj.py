# %%
import math
import numpy as np

RAD2DEG = 180/math.pi
DEG2RAD = math.pi/180
R_earth = 6371.0

# read src_rec_file (for abs, cs, cr)
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

def read_obj_file(filename):
    doc = open(filename,'r')
    doc_input = doc.readlines()
    doc.close()

    all_obj = []
    for iline,info in enumerate(doc_input):

        if(iline == 0):
            continue

        tmp=info.split(",")
        obj = []

        # from id 2 to 7
        for i in range(2,7):
            obj.append(float(tmp[i]))

        all_obj.append(obj)

    return all_obj

# %%
fname = 'src_rec_obs.dat'
time_true,ev_time_true = read_src_rec_file(fname)

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_abs'
time_no_swap_abs,ev_time_no_swap_abs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_abs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_abs_cr'
time_no_swap_abs_cr,ev_time_no_swap_abs_cr = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_abs_cr = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_abs_cs'
time_no_swap_abs_cs,ev_time_no_swap_abs_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_abs_cs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_abs_cr_cs'
time_no_swap_abs_cr_cs,ev_time_no_swap_abs_cr_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_abs_cr_cs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_no_data'
time_no_swap_no_data,ev_time_no_swap_no_data = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_no_data = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_cr'
time_no_swap_cr,ev_time_no_swap_cr = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_cr = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_cs'
time_no_swap_cs,ev_time_no_swap_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_cs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_no_swap_cr_cs'
time_no_swap_cr_cs,ev_time_no_swap_cr_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_no_swap_cr_cs = read_obj_file("%s/objective_function.txt"%(fname))


fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_abs'
time_swap_abs,ev_time_swap_abs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_abs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_abs_cr'
time_swap_abs_cr,ev_time_swap_abs_cr = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_abs_cr = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_abs_cs'
time_swap_abs_cs,ev_time_swap_abs_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_abs_cs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_abs_cr_cs'
time_swap_abs_cr_cs,ev_time_swap_abs_cr_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_abs_cr_cs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_no_data'
time_swap_no_data,ev_time_swap_no_data = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_no_data = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_cr'
time_swap_cr,ev_time_swap_cr = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_cr = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_cs'
time_swap_cs,ev_time_swap_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_cs = read_obj_file("%s/objective_function.txt"%(fname))

fname = 'OUTPUT_FILES/OUTPUT_FILES_swap_cr_cs'
time_swap_cr_cs,ev_time_swap_cr_cs = read_src_rec_file("%s/src_rec_file_step_0000.dat"%(fname))
obj_swap_cr_cs = read_obj_file("%s/objective_function.txt"%(fname))

# %%
# compare the src_rec output file

error_no_swap_serial = ((time_no_swap_abs_cr_cs-time_true))
if (error_no_swap_serial.std() < 0.1):
    print("For time_no_swap_abs_cr_cs,  src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))
else:
    print("For time_no_swap_abs_cr_cs,  output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_no_swap_parallel = ((time_no_swap_abs_cr-time_true))
if (error_no_swap_parallel.std() < 0.1):
    print("For time_no_swap_abs_cr,     src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_parallel.std())))
else:
    print("For time_no_swap_abs_cr,     output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))


error_swap_serial = ((time_no_swap_abs_cs-time_true))
if (error_swap_serial.std() < 0.1):
    print("For time_no_swap_abs_cs,     src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_serial.std())))
else:
    print("For time_no_swap_abs_cs,     output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_swap_parallel = ((time_no_swap_abs-time_true))
if (error_swap_parallel.std() < 0.1):
    print("For time_no_swap_abs,        src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_parallel.std())))
else:
    print("For time_no_swap_abs,        output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_no_swap_serial = ((time_no_swap_cr_cs-time_true))
if (error_no_swap_serial.std() < 0.1):
    print("For time_no_swap_cr_cs,      src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))
else:
    print("For time_no_swap_cr_cs,      output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_no_swap_parallel = ((time_no_swap_cr-time_true))
if (error_no_swap_parallel.std() < 0.1):
    print("For time_no_swap_cr,         src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_parallel.std())))
else:
    print("For time_no_swap_cr,         output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))


error_swap_serial = ((time_no_swap_cs-time_true))
if (error_swap_serial.std() < 0.1):
    print("For time_no_swap_cs,         src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_serial.std())))
else:
    print("For time_no_swap_cs,         output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_swap_parallel = ((time_no_swap_no_data-time_true))
if (error_swap_parallel.std() < 0.1):
    print("For time_no_swap_no_data,    src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_parallel.std())))
else:
    print("For time_no_swap_no_data,    output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))


###############################  swap cases #####################

error_no_swap_serial = ((time_swap_abs_cr_cs-time_true))
if (error_no_swap_serial.std() < 0.1):
    print("For time_swap_abs_cr_cs,     src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))
else:
    print("For time_swap_abs_cr_cs,     output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_no_swap_parallel = ((time_swap_abs_cr-time_true))
if (error_no_swap_parallel.std() < 0.1):
    print("For time_swap_abs_cr,        src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_parallel.std())))
else:
    print("For time_swap_abs_cr,        output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))


error_swap_serial = ((time_swap_abs_cs-time_true))
if (error_swap_serial.std() < 0.1):
    print("For time_swap_abs_cs,        src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_serial.std())))
else:
    print("For time_swap_abs_cs,        output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_swap_parallel = ((time_swap_abs-time_true))
if (error_swap_parallel.std() < 0.1):
    print("For time_swap_abs,           src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_parallel.std())))
else:
    print("For time_swap_abs,           output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_no_swap_serial = ((time_swap_cr_cs-time_true))
if (error_no_swap_serial.std() < 0.1):
    print("For time_swap_cr_cs,         src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))
else:
    print("For time_swap_cr_cs,         output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_no_swap_parallel = ((time_swap_cr-time_true))
if (error_no_swap_parallel.std() < 0.1):
    print("For time_swap_cr,            src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_no_swap_parallel.std())))
else:
    print("For time_swap_cr,            output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))


error_swap_serial = ((time_swap_cs-time_true))
if (error_swap_serial.std() < 0.1):
    print("For time_swap_cs,            src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_serial.std())))
else:
    print("For time_swap_cs,            output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))

error_swap_parallel = ((time_swap_no_data-time_true))
if (error_swap_parallel.std() < 0.1):
    print("For time_swap_no_data,       src_rec.dat read and output are CORRECT! (match coeffient: %s)"%(1-abs(error_swap_parallel.std())))
else:
    print("For time_swap_no_data,       output src_rec.dat is not CONSISTENT with the input src_rec.dat, please check the output src_rec file and the code! (match coeffient: %s)"%(1-abs(error_no_swap_serial.std())))


# %%
# check obj output

print("The output obj files are list here: ")
print("  Use data types        total obj         obj of abs             obj of cs             obj of cr")
tmp = obj_no_swap_abs_cr_cs
print("no_swap_abs_cr_cs: %8.2f(total) %8.2f(    use abs) %8.2f(    use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_abs_cr
print("no_swap_abs_cr:    %8.2f(total) %8.2f(    use abs) %8.2f(not use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_abs_cs
print("no_swap_abs_cs:    %8.2f(total) %8.2f(    use abs) %8.2f(    use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_abs
print("no_swap_abs:       %8.2f(total) %8.2f(    use abs) %8.2f(not use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_cr_cs
print("no_swap_cr_cs:     %8.2f(total) %8.2f(not use abs) %8.2f(    use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_cr
print("no_swap_cr:        %8.2f(total) %8.2f(not use abs) %8.2f(not use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_cs
print("no_swap_cs:        %8.2f(total) %8.2f(not use abs) %8.2f(    use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_no_data
print("no_swap_no_data:   %8.2f(total) %8.2f(not use abs) %8.2f(not use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
print("")
tmp = obj_no_swap_abs_cr_cs
print("swap_abs_cr_cs:    %8.2f(total) %8.2f(    use abs) %8.2f(    use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_abs_cr
print("swap_abs_cr:       %8.2f(total) %8.2f(    use abs) %8.2f(not use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_abs_cs
print("swap_abs_cs:       %8.2f(total) %8.2f(    use abs) %8.2f(    use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_abs
print("swap_abs:          %8.2f(total) %8.2f(    use abs) %8.2f(not use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_cr_cs
print("swap_cr_cs:        %8.2f(total) %8.2f(not use abs) %8.2f(    use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_cr
print("swap_cr:           %8.2f(total) %8.2f(not use abs) %8.2f(not use cs ) %8.2f(    use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_cs
print("swap_cs:           %8.2f(total) %8.2f(not use abs) %8.2f(    use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))
tmp = obj_no_swap_no_data
print("swap_no_data:      %8.2f(total) %8.2f(not use abs) %8.2f(not use cs ) %8.2f(not use cr )"%(tmp[0][0],tmp[0][1],tmp[0][2],tmp[0][3]))


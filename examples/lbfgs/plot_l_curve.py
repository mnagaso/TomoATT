import numpy as np
import matplotlib.pyplot as plt

list_obj_files = [
    "./OUTPUT_FILES/objective_function_l_0.00001.txt",
    "./OUTPUT_FILES/objective_function_l_0.0001.txt",
    "./OUTPUT_FILES/objective_function_l_0.001.txt",
    "./OUTPUT_FILES/objective_function_l_0.01.txt",
    "./OUTPUT_FILES/objective_function_l_0.1.txt",
    "./OUTPUT_FILES/objective_function_l_1.0.txt",
]

if __name__ == "__main__":

    list_fit = []  # fitting error
    list_regl = [] # regularization term
    list_coef = [] # regularization coefficient

    # find the last row which has 1 for the last column
    for file in list_obj_files:
        with open(file, 'r') as f:
            try:
                lines = f.readlines()
                obj_line= 2
                for i, line in enumerate(lines):
                    if line.strip().split(",")[-1] == "1":
                        break
                    else:
                        obj_line = i

                # split the line
                line = lines[obj_line].strip().split(",")

                # append the values
                lambda_val = float(file.split("_")[-1].rstrip(".txt"))
                obj = float(line[3])
                term_regl_tmp = float(line[5]) #/ obj
                term_regl = term_regl_tmp*2/lambda_val
                term_fit = (obj-term_regl_tmp) #/ obj
                list_fit.append(term_fit)
                list_regl.append(term_regl)
                list_coef.append(lambda_val)
            except:
                print("Error reading file: ", file)

    # plot the L-curve
    plt.figure()
    plt.loglog(list_fit, list_regl, 'o-')
    plt.xlabel("residual norm (||Ax_i-b||^2)")
    plt.ylabel("solution norm (||x_i-x_0||^2)")
    # add the regularization coefficient
    for i in range(len(list_coef)):
        if i == 0:
            plt.text(list_fit[i], list_regl[i], "lambda="+str(list_coef[i]))
        else:
            plt.text(list_fit[i], list_regl[i], str(list_coef[i]))
    # tick labels in scientific notation
    plt.grid()
    plt.show()

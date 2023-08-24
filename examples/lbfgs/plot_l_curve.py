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

    # read only the second row of the files
    for file in list_obj_files:
        with open(file, 'r') as f:
            lines = f.readlines()
            print(lines[2].strip())

            # split the line
            line = lines[2].strip().split(",")

            # append the values
            lambda_val = float(file.split("_")[-1].rstrip(".txt"))
            list_fit.append(float(line[3])-float(line[5]))
            list_regl.append(float(line[5]))
            list_coef.append(lambda_val)


    # plot the L-curve
    plt.figure()
    plt.loglog(list_fit, list_regl, 'o-')
    plt.xlabel("residual norm (||Ax_i-b||^2)")
    plt.ylabel("solution norm (0.5*lambda*||x_i-x_0||^2)")
    # add the regularization coefficient
    for i in range(len(list_coef)):
        if i == 0:
            plt.text(list_fit[i], list_regl[i], "lambda="+str(list_coef[i]))
        else:
            plt.text(list_fit[i], list_regl[i], str(list_coef[i]))
    # tick labels in scientific notation
    plt.grid()
    plt.show()

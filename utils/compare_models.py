#!/bin/python3

import os
import h5py
import numpy as np
import argparse

def read_model_h5(model_path):

    with h5py.File(model_path, 'r') as f:
        # read xi eta vel
        xi = f['xi'][:]
        eta = f['eta'][:]
        vel = f['vel'][:]

    return xi, eta, vel

# main function
if __name__ == '__main__':
    """
    Compare two models
    parth args
    -t : true model file path
    -r : result model file path
    """

    parser = argparse.ArgumentParser(description='Compare two models')
    parser.add_argument('-t', '--true', type=str, help='true model file path')
    parser.add_argument('-r', '--result', type=str, help='result model file path')
    args = parser.parse_args()

    # read true model
    xi_true, eta_true, vel_true = read_model_h5(args.true)
    # read result model
    xi_result, eta_result, vel_result = read_model_h5(args.result)

    # check shapes
    if xi_true.shape != xi_result.shape:
        print('xi shape not match')
        exit(1)

    if eta_true.shape != eta_result.shape:
        print('eta shape not match')
        exit(1)

    if vel_true.shape != vel_result.shape:
        print('vel shape not match')
        exit(1)

    # print model info
    print('model info: ')
    print('vel shape: ', vel_true.shape)

    # compare
    print('vel max error: ', np.max(np.abs(vel_true - vel_result)))
    print('xi max error: ', np.max(np.abs(xi_true - xi_result)))
    print('eta max error: ', np.max(np.abs(eta_true - eta_result)))

    # L2 norm
    print('vel L2 norm: ', np.linalg.norm(vel_true - vel_result))
    print('xi L2 norm: ', np.linalg.norm(xi_true - xi_result))
    print('eta L2 norm: ', np.linalg.norm(eta_true - eta_result))

    # exit
    exit(0)



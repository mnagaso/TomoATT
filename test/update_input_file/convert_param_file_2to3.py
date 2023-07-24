import argparse
from ruamel.yaml import YAML
from contextlib import suppress

def map_value_to_bool(params_in, target_key, orig_key=None):
    """
    Map an integer value to a boolean and update the target key in params_in.

    Parameters:
        params_in (dict): The input dictionary.
        target_key (str): The key whose value needs to be mapped to a boolean.

    Returns:
        None
    """
    if orig_key is None:
        orig_key = target_key
        print('key name {} is changed to {}'.format(orig_key, target_key))

    value = params_in.get(orig_key, None)
    if value is not None:
        params_in[target_key] = bool(value)
        print('value {} type is changed from int to bool'.format(target_key))

    # remove the old key
    if orig_key != target_key:
        params_in.pop(orig_key, None)


def move_value(params_out, params_in):
    try:
        params_out = params_in
    except KeyError:
        print('Cannot move value from {} to {}'.format(params_in, params_out))


if __name__ == '__main__':
    # parse the argument for the input file
    parser = argparse.ArgumentParser(description='Convert a parameter file from version 2 to version 3')
    parser.add_argument('-i', '--input', help='Input file name', required=True)

    args = parser.parse_args()
    infile = args.input

    # path to the v3 model file
    v3_model_file = './params_model_v3.yaml'

    yaml = YAML()

    # read the input file
    try:
        with open(infile, 'r') as f:
            str_in = f.read()
            params_in = yaml.load(str_in)

    except IOError:
        raise ValueError('Cannot read the input file')

    # read the v3 model file
    try:
        with open(v3_model_file, 'r') as f:
            str_v3 = f.read()
            params_v3 = yaml.load(str_v3)
    except IOError:
        raise ValueError('Cannot read the v3 model file')

    # check the version of the input file
    if params_in['version'] != 2:
        raise ValueError('The input file is not version 2')

    # change version to 3
    params_v3['version'] = 3

    # copy the values in the input file to the output file
    #
    # domain section
    #
    params_v3['domain'] = params_in['domain']

    #
    # source section
    #
    params_v3['source'] = params_in['source']

    #
    # model section
    #
    params_v3['model'] = params_in['model']

    #
    # parallel section
    #
    params_v3['parallel'] = params_in['parallel']

    # change parallel->use_gpu from 0,1 to false,true
    map_value_to_bool(params_v3['parallel'], 'use_gpu')

    #
    # output_setting section
    #
    move_value(params_v3['output_setting']['output_dir'],params_in['inversion']['output_dir'])

    map_value_to_bool(params_v3['output_setting'], 'output_source_field', 'is_output_source_field')
    map_value_to_bool(params_v3['output_setting'], 'output_model_dat', 'is_output_model_dat')
    map_value_to_bool(params_v3['output_setting'], 'output_final_model', 'is_output_final_model')
    map_value_to_bool(params_v3['output_setting'], 'output_in_process', 'is_output_in_process')
    map_value_to_bool(params_v3['output_setting'], 'single_precision_output', 'is_single_precision_output')

    # remove the old key 'output_setting'->'is_verbose_output'
    params_v3['output_setting'].pop('is_verbose_output', None)

    move_value(params_v3['output_setting']['output_file_format'], params_in['calculation']['output_file_format'])

    #
    # run_mode section
    #
    move_value(params_v3['run_mode'], params_in['inversion']['run_mode'])

    #
    # model_update section
    #
    move_value(params_v3['model_update']['max_iterations'],                            params_in['inversion']['max_iterations_inv'])
    move_value(params_v3['model_update']['optim_method'],                              params_in['inversion']['optim_method'])
    move_value(params_v3['model_update']['step_length'],                               params_in['inversion']['step_size'])
    move_value(params_v3['model_update']['optim_method_0']['step_length_decay'],       params_in['inversion']['step_size_decay'])
    move_value(params_v3['model_update']['optim_method_0']['step_length_sc'],          params_in['inversion']['step_size_sc'])
    move_value(params_v3['model_update']['optim_method_1_2']['max_sub_iterations'],    params_in['inversion']['max_sub_iterations'])
    move_value(params_v3['model_update']['optim_method_1_2']['regularization_weight'], params_in['inversion']['regularization_weight'])
    move_value(params_v3['model_update']['smoothing']['smooth_method'],                params_in['inversion']['smooth_method'])
    move_value(params_v3['model_update']['smoothing']['l_smooth_rtp'],                 params_in['inversion']['l_smooth_rtp'])
    move_value(params_v3['model_update']['n_inversion_grid'],                          params_in['inversion']['n_inversion_grid'])
    with suppress(KeyError) : move_value(params_v3['model_update']['type_invgrid_dep'],                          params_in['inversion']['type_dep_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['type_invgrid_lat'],                          params_in['inversion']['type_lat_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['type_invgrid_lon'],                          params_in['inversion']['type_lon_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['n_inv_dep_lat_lon'],                         params_in['inversion']['n_inv_dep_lat_lon'])
    with suppress(KeyError) : move_value(params_v3['model_update']['min_max_dep_inv'],                           params_in['inversion']['min_max_dep_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['min_max_lat_inv'],                           params_in['inversion']['min_max_lat_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['min_max_lon_inv'],                           params_in['inversion']['min_max_lon_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['dep_inv'],                                   params_in['inversion']['dep_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['lat_inv'],                                   params_in['inversion']['lat_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['lon_inv'],                                   params_in['inversion']['lon_inv'])
    with suppress(KeyError) : move_value(params_v3['model_update']['sta_correction_file'],                       params_in['inversion']['sta_correction_file'])
    with suppress(KeyError) : move_value(params_v3['model_update']['update_slowness'],                          params_in['inversion']['is_inv_slowness'])
    with suppress(KeyError) : move_value(params_v3['model_update']['update_azi_ani'],                           params_in['inversion']['is_inv_azi_ani'])
    with suppress(KeyError) : move_value(params_v3['model_update']['update_rad_ani'],                           params_in['inversion']['is_inv_rad_ani'])
    map_value_to_bool(params_v3['model_update'], 'update_slowness')
    map_value_to_bool(params_v3['model_update'], 'update_azi_ani')
    map_value_to_bool(params_v3['model_update'], 'update_rad_ani')

    with suppress(KeyError) : move_value(params_v3['model_update']['depth_taper'],                              params_in['inversion']['kernel_taper'])

    #
    # relocation section
    #
    # replocation section is completely new in v3, so we don't need to move any value

    #
    # inversion_strategy section
    #
    # inversion_strategy section is completely new in v3, so we don't need to move any value

    #
    # calculation section
    #
    move_value(params_v3['calculation'], params_in['calculation'])


    # write the output file with adding .v3 to the file name
    outfile = infile + '.v3.yaml'
    with open(outfile, 'w') as f:
        yaml.dump(params_v3, f)







import argparse
from ruamel.yaml import YAML


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
        print ('value {} type is changed from int to bool'.format(target_key))

    # remove the old key
    if orig_key != target_key:
        params_in.pop(orig_key, None)

def move_value(params_in, params_out):
    try:
        params_out = params_in
    except:
        print('Cannot move value from {} to {}'.format(params_in, params_out))

if __name__ == '__main__':
    # parse the argument for the input file
    parser = argparse.ArgumentParser(description='Convert a parameter file from version 2 to version 3')
    parser.add_argument('-i', '--input', help='Input file name', required=True)

    args = parser.parse_args()
    infile = args.input

    # path to the v3 model file
    v3_model_file = 'test/update_input_file/params_model_v3.yaml'

    yaml = YAML()

    # read the input file
    try:
        with open(infile, 'r') as f:
            str_in = f.read()
            params_in = yaml.load(str_in)

    except:
        raise ValueError('Cannot read the input file')

    # read the v3 model file
    try:
        with open(v3_model_file, 'r') as f:
            str_v3 = f.read()
            params_v3 = yaml.load(str_v3)
    except:
        raise ValueError('Cannot read the v3 model file')

    # check the version of the input file
    if params_in['version'] != 2:
        raise ValueError('The input file is not version 2')

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
    move_value(params_v3['run_mode'], params_in['inversoin']['run_mode'])

    #
    # model_update section
    #

    # write the output file with adding .v3 to the file name
    outfile = infile + '.v3.yaml'
    with open(outfile, 'w') as f:
        yaml.dump(params_v3, f)







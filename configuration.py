import ConfigParser

__author__ = 'quepas, jakubsieradzki'

def build_default():
    return {
                # Local mean (1), expectation-maximization (2)
                'ex_filter_type': 1,
                'ex_window_size': 3,
                # number of iterations of the EM algorithm
                'ex_iterations': 10,
                # sigma for LPF filter
                'lpf_f': 3.4,
                # sigma for LPF filter used to smooth sigma(x) in SNR
                'lpf_f_SNR': 1.2,
                # sigma for LPF filter used to smooth Rician corrected noise map
                'lpf_f_Rice': 5.4,
                'input_filename': 'data/input/noisy.csv',
                'input_filename_SNR': 'data/input/snr.csv',
                'output_filename_Gaussian': 'data/output/gaussian_map.csv',
                'output_filename_Rician': 'data/output/rician_map.csv'
            }

def build_from_file(filename):
    config_parser = ConfigParser.ConfigParser()
    config_parser.read(filename)
    config = build_default()
    config['ex_filter_type'] = config_parser.getint("config", "ex_filter_type")
    config['ex_window_size'] = config_parser.getint("config", "ex_window_size")
    config['ex_iterations'] = config_parser.getint("config", "ex_iterations")
    config['lpf_f'] = config_parser.getfloat("config", "lpf_f")
    config['lpf_f_SNR'] = config_parser.getfloat("config", "lpf_f_SNR")
    config['lpf_f_Rice'] = config_parser.getfloat("config", "lpf_f_Rice")
    config['input_filename'] = config_parser.get("config", "input_filename")
    if config_parser.has_option("config", "input_filename_SNR"):
        config['input_filename_SNR'] = config_parser.get("config", "input_filename_SNR")
    else:
        config['input_filename_SNR'] = ""
    config['output_filename_Gaussian'] = config_parser.get("config", "output_filename_Gaussian")
    config['output_filename_Rician'] = config_parser.get("config", "output_filename_Rician")
    return config

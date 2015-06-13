from numpy import genfromtxt, array, savetxt
from time import clock
from algorithm import rice_homomorf_est

__author__ = 'quepas, jakubsieradzki'

def run_single_experiment(config):
    now = clock()
    noisy = genfromtxt(config['input_filename'], delimiter=',')
    filename_snr = config['input_filename_SNR']
    if not filename_snr:
        snr = array([0])
    else:
        snr = genfromtxt(filename_snr, delimiter=',')
    (rician_map, gaussian_map) = rice_homomorf_est(noisy, snr, config['lpf_f'], config['ex_filter_type'], config=config);
    savetxt(config['output_filename_Gaussian'], gaussian_map, delimiter=',')
    savetxt(config['output_filename_Rician'], rician_map, delimiter=',')
    print "Run time : " + str(clock() - now) + "s"

def run_example_experiments(config):
    input_noisy = genfromtxt(config['input_filename'], delimiter=',')
    input_snr = genfromtxt(config['input_filename_SNR'], delimiter=',')
    lpf = config['lpf_f']

    # Expectation-maximization with known SNR
    now = clock()
    (rice_map, gaussian_map) = rice_homomorf_est(input_noisy, input_snr, lpf, 2);
    savetxt("data/output/rician_map_em_snr.csv", rice_map, delimiter=',')
    savetxt("data/output/gaussian_map_em_snr.csv", gaussian_map, delimiter=',')
    print "Expectation-maximization with known SNR time: " + str(clock() - now) + "s"

    # Expectation-maximization with unknown SNR
    now = clock()
    (rice_map, gaussian_map) = rice_homomorf_est(input_noisy, array([0]), lpf, 2);
    savetxt("data/output/rician_map_em.csv", rice_map, delimiter=',')
    savetxt("data/output/gaussian_map_em.csv", gaussian_map, delimiter=',')
    print "Expectation-maximization with unknown SNR time: " + str(clock() - now) + "s"

    # Local mean estimation with known SNR
    now = clock()
    (rice_map, gaussian_map) = rice_homomorf_est(input_noisy, input_snr, lpf, 1);
    savetxt("data/output/rician_map_lm_snr.csv", rice_map, delimiter=',')
    savetxt("data/output/gaussian_map_lm_snr.csv", gaussian_map, delimiter=',')
    print "Local mean estimation with known SNR time: " + str(clock() - now) + "s"

    # Local mean estimation with unknown SNR
    now = clock()
    (rice_map, gaussian_map) = rice_homomorf_est(input_noisy, array([0]), lpf, 1);
    savetxt("data/output/rician_map_lm.csv", rice_map, delimiter=',')
    savetxt("data/output/gaussian_map_lm.csv", gaussian_map, delimiter=',')
    print "Local mean estimation with unknown SNR time: " + str(clock() - now) + "s"

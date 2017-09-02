'''
filtering lib
'''
import numbers
import numpy as np
import utils as utl
import matplotlib.pyplot as plt

'''
Python implementation of lanczos filtering
'''
def lanczos_filter(time_series, spl_intvl=1, cut_off_freq=None,\
        num_of_coeffs=100, mode='low'):

    ''' input checking '''

    # time_series has to be a one dimensional numpy array
    assert (isinstance(time_series, (np.ndarray)) & (time_series.ndim == 1)),\
            'Error: Incorrect input dimension.'

    # time_series has to be a real vector
    assert np.isrealobj(time_series),\
            'Error: The input time series should be a real object.'

    # spl_intvl (sampling interval) has to be a real number
    assert isinstance(spl_intvl, numbers.Real),\
            'Error: The input sampling interval has to be a real number.'
    nyquist_freq = 1.0 / (2 * spl_intvl)

    # cut off freq has to be a positive real number
    if cut_off_freq is None :
        print("Half Nyquist Frequency is used as cut off frequency.\n")
        cut_off_freq = nyquist_freq / 2.0
    assert isinstance(cut_off_freq, numbers.Real) & (cut_off_freq > 0),\
            'Error: Cut off frequency has to be a positive number.'

    # number of coefficients has to be a positive integer
    assert isinstance(num_of_coeffs, numbers.Integral) & (num_of_coeffs > 0),\
            'Error: number of coefficients has to be a positive integer.'

    # mode has to be either 'high' or 'low'
    assert isinstance(mode, basestring) & ((mode == 'low') | (mode == 'high')),\
            'Error: Only high/low pass filter are supported.'

    if mode == 'high':
        LoH = 1
    else :
        LoH = 0

    ''' Start '''

    # Normalize the cut off frequency with the Nyquist frequency
    cut_off_freq = cut_off_freq / nyquist_freq
    lanczos_coeff = utl.lanczos_filter_coef(cut_off_freq, num_of_coeffs)
    lanczos_coeff = lanczos_coeff[:,LoH]

    # get filter in frequency space
    window, freq = utl.spectral_window(lanczos_coeff, len(time_series))
    freq = freq * nyquist_freq

    # Replace NaN's with the mean (ideas?):
    index_NaN = np.isnan(time_series)
    ts_mean = np.mean(time_series[~index_NaN])
    time_series[index_NaN] = ts_mean

    # Filtering
    filtered_ts, fft_ts = utl.spectral_filtering(time_series, window)

    return filtered_ts, lanczos_coeff, window, fft_ts, freq

    '''
    Python implementation of eval_KZ_choice
    '''
def eval_kz_choice(length, kz_choice, mode='low'):

    ''' input checking '''
    assert isinstance(kz_choice,(np.ndarray)) & (kz_choice.shape == (2,)), \
            'Error: kz_choice has to be a two element 1d array'

    ''' create a delta function '''
    delta = np.zeros(length)

    if (length % 2 == 0) : #even
        delta[int(length / 2.0)] = 1e8
    else : # odd
        delta[np.ceil(length / 2.) - 1] = 1e8

    freq, power_delta = utl.power_spectral(delta, len(delta))

    mov_avg, t = utl.kz_low_pass(delta, kz_choice[0], kz_choice[1])

    if mode == 'high':
        mov_avg = delta[t] - mov_avg

    freq2, power_mov = utl.power_spectral(mov_avg, len(delta))

    cutperiod = length / (np.argmin(np.absolute(power_delta[1] * np.exp(-1) - \
        power_mov)) + 1)

    plt.show()
    return cutperiod









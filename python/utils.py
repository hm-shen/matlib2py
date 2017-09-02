import numbers
import numpy as np

def lowpass_cosine_filter_coeff(cut_off_freq, num_of_coeffs):
    tmp = np.linspace(1, num_of_coeffs, num_of_coeffs)
    rst = np.concatenate((np.array([1.0]), np.sin(np.pi * tmp * cut_off_freq) \
                / (np.pi * tmp * cut_off_freq))) * cut_off_freq
    return rst

def cmpt_sigma_factors(cut_off_freq, num_of_coeffs):
    tmp = np.linspace(1, num_of_coeffs, num_of_coeffs)
    rst = np.concatenate((np.array([1.0]), np.sin(np.pi * tmp / num_of_coeffs) \
            / (np.pi * tmp / num_of_coeffs)))
    return rst

def lanczos_filter_coef(cut_off_freq, num_of_coeffs):

    ''' compute positive coefficients of Lanczos [low high] pass '''
    # compute lowpass coeffs
    lowpass_coeff = lowpass_cosine_filter_coeff(cut_off_freq, \
            num_of_coeffs)
    sigma = cmpt_sigma_factors(cut_off_freq, num_of_coeffs)
    hk_low = lowpass_coeff * sigma

    # compute highpass coeffs
    hk_high = -1 * hk_low
    hk_high[0] = hk_high[0] + 1

    # concatenate hk_low and hk_high
    return np.concatenate((hk_low[:,np.newaxis], hk_high[:,np.newaxis]),axis=1)

def spectral_window(coeffs, len_of_ts):

    # compute window in frequency space
    freq = np.arange(0, 1+1e-15, (2.0 / len_of_ts))
    window = np.zeros(freq.shape)
    tmp = np.linspace(1, len(coeffs)-1, len(coeffs)-1)
    for index in range(len(freq)):
        window[index] = coeffs[0] + 2 * \
                np.sum(coeffs[1 :] * np.cos(tmp * np.pi * freq[index]))
    return window, freq

def spectral_filtering(time_series, window):
    # Filtering in frequency space is multiplication,
    # (convolution in time space).

    series_length = len(time_series)
    fft_ts = np.fft.fft(time_series)
    fft_ts = fft_ts[0 : int(np.floor(series_length / 2) + 1)]
    fft_ts_filtered = fft_ts * window
    fft_ts_filtered = np.concatenate((fft_ts_filtered, np.conjugate(\
            fft_ts_filtered[series_length - len(fft_ts_filtered) : 0 : -1])))
    ts_filtered = np.real(np.fft.ifft(fft_ts_filtered))

    return ts_filtered, fft_ts


import numbers
import numpy as np
import matplotlib.pyplot as plt


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

def moving_average(time_series, num_of_points):

    ''' init parameters '''
    begin_ind = int(np.ceil(num_of_points / 2.)) - 1
    end_ind = len(time_series)- int(np.floor(num_of_points / 2.) + 1)
    selected_ind = np.linspace(begin_ind, end_ind, (end_ind - begin_ind)+1)

    ''' calculate moving sum using cumsum '''
    cumsum_of_ts = np.cumsum(np.insert(time_series,0,0))
    move_avg = (cumsum_of_ts[num_of_points:] -\
            cumsum_of_ts[:-num_of_points]) / float(num_of_points)

    return move_avg, selected_ind

def kz_low_pass(data, m_points, k_iters):

    t = np.linspace(1,len(data), len(data))
    mov_avg = data

    for index in range(k_iters):
        inner_mov_avg, tmp = moving_average(mov_avg, m_points)
        t = t[0] + tmp
        mov_avg = inner_mov_avg

    return mov_avg, (t - 1).astype(int)

def power_spectral(signal, num_of_points, plot_flag=True):

    # fft
    sig_freq = np.fft.fft(signal, num_of_points)
    power_sig = sig_freq * np.conjugate(sig_freq) / num_of_points
    freq = len(signal) * np.linspace(0,num_of_points/2,1+num_of_points/2) /  num_of_points

    if plot_flag:
        plt.figure()
        plt.plot(freq, power_sig[: int(np.floor(num_of_points/2.) + 1)])
        plt.title('Frequency content of signal')
        plt.xlabel('frequency (Hz)')

    return freq, power_sig


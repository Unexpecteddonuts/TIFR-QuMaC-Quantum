import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import cauchy
from scipy.signal import find_peaks
import json, os

def cauchy_fit(x, x0, gamma, A, y0):
    return y0 + A * cauchy.pdf(x - x0, scale=gamma)

def ext_BW(freq, ydata, noise_threshold, structure_threshold, min_averages, max_averages):
    # Find peaks
    peaks, _ = find_peaks(ydata)
    peak_indexes = peaks[np.argsort(peaks[ydata[peaks]])[::-1]]
    xdata = freq[peak_indexes]
    ydata = ydata[peak_indexes]

    num_averages = detect_snr(ydata, noise_threshold, structure_threshold, min_averages, max_averages)

    # Calculate average and standard deviation of signal
    avg_signal = np.mean(ydata[:num_averages])
    std_signal = np.std(ydata[:num_averages])

    # Calculate SNR
    if std_signal > 0:
        snr = avg_signal / std_signal
    else:
        snr = 0

    # Initialize properties with None
    cavity_frequency = bandwidth = internal_bandwidth = external_bandwidth = internal_quality_factor = None

    if snr > noise_threshold and snr > structure_threshold:
        # Initial guess
        x0_guess = xdata[np.argmax(ydata)]
        y0_guess = np.min(ydata)
        A_guess = np.max(ydata) - y0_guess
        gamma_guess = (xdata[-1] - xdata[0]) / 2

        # Curve fit
        p0 = [x0_guess, gamma_guess, A_guess, y0_guess]
        params, cov = curve_fit(cauchy_fit, freq, ydata, p0=p0)
        x0, gamma, A, y0 = params

        # Calculate FWHM
        fwhm = 2 * gamma

        # Resonator properties
        bandwidth = fwhm
        cavity_frequency = x0
        internal_bandwidth = fwhm / (1 + A)
        external_bandwidth = fwhm - internal_bandwidth
        internal_quality_factor = cavity_frequency / internal_bandwidth

    properties = {
        'cavity_frequency': cavity_frequency,
        'bandwidth': bandwidth,
        'internal_bandwidth': internal_bandwidth,
        'external_bandwidth': external_bandwidth,
        'internal_quality_factor': internal_quality_factor,
        'num_averages': num_averages
    }

    return properties

def detect_snr(signal, noise_threshold, structure_threshold, min_averages, max_averages):
    num_averages = 0
    snr = 0

    while num_averages < max_averages:
        num_averages += 1

        # Calculate average and standard deviation of signal
        avg_signal = np.mean(signal[:num_averages])
        std_signal = np.std(signal[:num_averages])

        # Calculate SNR
        if std_signal > 0:
            snr = avg_signal / std_signal

        # Check if signal is present and has clear structure
        if snr > noise_threshold and snr > structure_threshold:
            # Check if minimum number of averages have been reached
            if num_averages >= min_averages:
                break

    return num_averages

# Defining power range
pmin, pmax = -100, 10
del_p = 10
p_vals = np.arange(pmin, pmax + del_p - 1, del_p)

# SNR and averaging parameters
noise_threshold = 3
structure_threshold = 5
min_averages = 10
max_averages = 100

# storage lists
kint_list, kext_list, Qint_list, Qext_list, fit_param_list, Q_with_pow_list, num_avgs_list = [], [], [], [], [], [], []

# Looping over power values and filepaths
for p in p_vals:
    # Creating the filepath for each power value
    filepath = f'F:\\Expt Data\\date\\folder\\data_{p}.0_dBm.txt'

    # Loading and transposing the data
    data = np.transpose(np.loadtxt(filepath))
    freq = data[0]
    ydata = data[1]

    # Calculating the fit parameters and properties
    properties = ext_BW(freq, ydata, noise_threshold, structure_threshold, min_averages, max_averages)

    fit_param_list.append([p, properties['cavity_frequency'], properties['bandwidth'], properties['internal_bandwidth'],
                           properties['external_bandwidth'], properties['internal_quality_factor']])
    
    num_avgs_list.append([p, properties['num_averages']])

    Q_with_pow_list.append([p, properties['cavity_frequency'], properties['internal_quality_factor'],
                            properties['internal_quality_factor_error'], properties['external_quality_factor'],
                            properties['external_quality_factor_error']])

    # Appending to list
    kint_list.append(properties['internal_bandwidth'])
    kext_list.append(properties['external_bandwidth'])
    Qint_list.append(properties['internal_quality_factor'])

# Saving the data to text files
np.savetxt(os.path.join(os.getcwd(), "Lorentian_fit_params.txt"), fit_param_list,
           header="Power Freq Bandwidth Internal_BW External_BW Internal_Qfactor")
np.savetxt(os.path.join(os.getcwd(), "Q_factors.txt"), Q_with_pow_list,
           header="Power Freq Internal_Qfactor Internal_QFactor_err External_QFactor External_QFactor_err")
np.savetxt(os.path.join(os.getcwd(), "Num_averages.txt"), num_avgs_list,
           header="Power Num_Averages")
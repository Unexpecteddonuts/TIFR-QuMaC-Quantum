import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import cauchy
from scipy.signal import find_peaks
import os

def cauchy_fit(x, x0, gamma, A, y0):
    return y0 + A * cauchy.pdf(x - x0, scale=gamma)

def ext_BW(freq, ydata):
    # Find peaks
    peaks, _ = find_peaks(ydata)
    peak_indexes = peaks[np.argsort(peaks[ydata[peaks]])[::-1]]
    xdata = freq[peak_indexes]
    ydata = ydata[peak_indexes]

    # Initial guess
    x0_guess = xdata[np.argmax(ydata)]
    y0_guess = np.min(ydata)
    A_guess = np.max(ydata) - y0_guess
    gamma_guess = (xdata[-1] - xdata[0]) / 2

    # Curve fit2
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
    external_quality_factor = cavity_frequency / external_bandwidth

    # Calculate errors
    fwhm_err = np.sqrt(cov[1, 1])
    cavity_frequency_err = np.sqrt(cov[0, 0])
    internal_bandwidth_err = fwhm_err / (1 + A)
    external_bandwidth_err = fwhm_err - internal_bandwidth_err
    internal_quality_factor_err = cavity_frequency_err / internal_bandwidth
    external_quality_factor_err = cavity_frequency_err / external_bandwidth
    
    
    properties = {
        'cavity_frequency': cavity_frequency,
        'bandwidth': bandwidth,
        'internal_bandwidth': internal_bandwidth,
        'external_bandwidth': external_bandwidth,
        'internal_quality_factor': internal_quality_factor,
        'internal_quality_factor_error': internal_quality_factor_err,
        'external_quality_factor': external_quality_factor,
        'external_quality_factor_error': external_quality_factor_err
    }

    return properties

# Defining power range
# pmin, pmax = -100, 10
# del_p = 10
# p_vals = np.arange(pmin,pmax + del_p -1, del_p)

# storage lists
kint_list, kext_list, Qint_list, Qext_list, fit_param_list, Q_with_pow_list = [], [], [], [], [], []

# Looping over power values and filepaths
# for p in p_vals:
# Creating the filepath for each power value
filepath = 'CavityResponse_forAakif.txt'

# Loading and transposing the data
data = np.transpose(np.loadtxt(filepath))
freq = data[0]
ydata = data[1]

# Calculating the fit parameters and properties
properties = ext_BW(freq, ydata)

fit_param_list.append([properties['cavity_frequency'], properties['bandwidth'], properties['internal_bandwidth'],
                        properties['external_bandwidth'], properties['internal_quality_factor']])

Q_with_pow_list.append([properties['cavity_frequency'], properties['internal_quality_factor'],
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
        header="Power Freq Internal_Qfactor Internal_QFactor_err External_QFactor External_QFactor_err ")
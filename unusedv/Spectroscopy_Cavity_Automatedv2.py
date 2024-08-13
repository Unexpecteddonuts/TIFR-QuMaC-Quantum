import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

def lorentzian(x,c,gam,a,y0):
    return y0 + (2*a/np.pi)*gam/(4*(x-c)**2 + gam**2)

def ext_BW(freq, ydata):
    
    # initial guess
    ymin = np.min(ydata)
    ymax = np.max(ydata)
    left_index = np.where(ydata > 0.5*(ymax+ymin))[0][0]    #FWHM
    right_index = np.where(ydata > 0.5*(ymax+ymin))[0][-1]   #FWHM
    BW_g = freq[right_index] - freq[left_index]
    wc_g = freq[np.argmax(ydata)]
    a_g = 0.5*np.pi*(ymax-ymin)*BW_g
    y0_g = ymin

    # curve fit
    res, cov = curve_fit(lorentzian,freq,ydata,[wc_g, BW_g, a_g, y0_g])
    
    f0, bw = res[0], res[1]
    f0_err, bw_err = np.sqrt(cov[0,0]), np.sqrt(cov[1,1])
    a, y0 = res[2], res[3]
    a_err, y0_err = np.sqrt(cov[2,2]), np.sqrt(cov[3,3])

    # calculation of resonator parameters
    H = 2*a/(np.pi*bw)
    H_err = H*(a_err/a + bw_err/bw)
    y_xc = H + y0
    y_xc_err = H_err + y0_err
    r = (y0 - y_xc)/(y0 + y_xc)
    r_err = r*(y0_err + y_xc_err)*(1/(-y0 + y_xc) + 1/(-y0 - y_xc))
    kint = bw/(1+r)
    kext = bw*r/(1+r)
    kint_err = kint*(bw_err/bw + r_err/(r+1))
    kext_err = bw_err + kint_err

    Qint = f0/kint
    Qext = f0/kext
    
    # Calculating errors
    Qint_err = Qint*(f0_err/f0 + kint_err/kint)
    Qext_err = Qext*(f0_err/f0 + kext_err/kext)

    f0 = np.round(f0*1e-9,10)
    f0_err = np.round(f0_err*1e-9,10)
    bw = np.round(bw*1e-9,10)
    bw_err = np.round(bw_err*1e-9,10)
    kint = np.round(kint*1e-9,10)
    kext = np.round(kext*1e-9,10)
    kint_err = np.round(kint_err*1e-9,10)
    kext_err = np.round(kext_err*1e-9,10)

    properties = {    #writing results to a dictionary
        'cavity_frequency': f0,
        'bandwidth': bw,
        'internal_bandwidth': kint,
        'external_bandwidth': kext,
        'internal_quality_factor': Qint,
        'internal_quality_factor_error': Qint_err,
        'external_quality_factor': Qext,
        'external_quality_factor_error': Qext_err
        }
            
            
    # Convert units and round to 1 kHz
    f0_ghz = round(properties['cavity_frequency'] * 1e-9, 3)
    bw_mhz = round(properties['bandwidth'] * 1e-6, 3)
    kint_mhz = round(properties['internal_bandwidth'] * 1e-6, 3)
    kext_mhz = round(properties['external_bandwidth'] * 1e-6, 3)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(freq, ydata, 'bo', label='Data')
    ax.plot(freq, lorentzian(freq, *res), 'r', label='Lorentzian Fit')
    ax.text(0.6, 0.9, 'Qint: {:.2f}'.format(Qint), horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Mearures Re(S11)')
    ax.set_title("plot")
    ax.legend()

    # Print results on the plot
    ax.text(0.05, 0.8, 'Cavity frequency: {:.3f} GHz'.format(f0_ghz), transform=ax.transAxes)
    ax.text(0.05, 0.75, 'Total bandwidth: {} MHz'.format(bw_mhz), transform=ax.transAxes)
    ax.text(0.05, 0.7, 'Internal bandwidth: {} MHz'.format(kint_mhz), transform=ax.transAxes)
    ax.text(0.05, 0.65, 'External bandwidth: {} MHz'.format(kext_mhz), transform=ax.transAxes)

    plt.show()
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
filepath = '2023-05-24_Real_0dBm.txt'

# Loading and transposing the data
data= np.transpose(np.loadtxt(filepath))
freq=data[0]
ydata=data[1]
print(freq, ydata)
# Calculating the fit parameters and properties
properties = ext_BW(freq, ydata)

# fit_param_list.append([p, properties['cavity_frequency'], properties['bandwidth'], properties['internal_bandwidth'],
#                         properties['external_bandwidth'], properties['internal_quality_factor']])

# Q_with_pow_list.append([p, properties['cavity_frequency'], properties['internal_quality_factor'], properties['internal_quality_factor_error'], properties['external_quality_factor'], properties['external_quality_factor_error']])

fit_param_list.append([properties['cavity_frequency'], properties['bandwidth'], properties['internal_bandwidth'],
                        properties['external_bandwidth'], properties['internal_quality_factor']])

Q_with_pow_list.append([properties['cavity_frequency'], properties['internal_quality_factor'], properties['internal_quality_factor_error'], properties['external_quality_factor'], properties['external_quality_factor_error']])

# Appending to list
kint_list.append(properties['internal_bandwidth'])
kext_list.append(properties['external_bandwidth'])

Qint_list.append(properties['internal_quality_factor'])

# Saving the data to text files
# np.savetxt(os.path.join(os.getcwd(), "Lorentian_fit_params.txt"), fit_param_list, 
#            header="Power Freq Bandwidth Internal_BW External_BW Internal_Qfactor")
# np.savetxt(os.path.join(os.getcwd(), "Q_factors.txt"), Q_with_pow_list, 
#            header="Power Freq Internal_Qfactor Internal_QFactor_err External_QFactor External_QFactor_err ")

np.savetxt(os.path.join(os.getcwd(), "Lorentian_fit_params.txt"), fit_param_list, 
           header="Freq Bandwidth Internal_BW External_BW Internal_Qfactor")
np.savetxt(os.path.join(os.getcwd(), "Q_factors.txt"), Q_with_pow_list, 
           header="Freq Internal_Qfactor Internal_QFactor_err External_QFactor External_QFactor_err ")
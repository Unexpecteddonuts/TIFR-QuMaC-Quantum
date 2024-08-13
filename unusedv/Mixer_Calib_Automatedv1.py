from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv2 import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
import os


simulate = False 
check_e_delay = False 

rr_no = 4
rr = f"rr{rr_no}"
out = adc_mapping[rr]
ro_len = ro_len_clk[str(rr_no)]
rep_rate_clk = 2500
rr_LO = config["elements"][rr]["mixInputs"]["lo_frequency"]

f_min = 10e6
f_max = 70e6
initial_df = 0.05e6
df_inc_factor = 1.2


if check_e_delay:
    f_min = -250e6
    f_max = 250e6
    initial_df = 1e6

initial_f_min = f_min
initial_f_max = f_max
qmm = QuantumMachinesManager()

# Finding FWHM and new range loop
for i in range(5):
    print("Running experiment")
    
    with program() as rr_spec:
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)

        with for_(n, 0, n < 200, n + 1):
            with for_(f, f_min, f < f_max, f + initial_df / df_inc_factor ** i):

                update_frequency(rr, f)
                wait(rep_rate_clk - ro_len, rr)
                measure("readout", rr, None,
                        demod.full("integW_cos", I, out),
                        demod.full("integW_minus_sin", Q, out))

                save(I, I_st)
                save(Q, Q_st)

        with stream_processing():
            I_st.buffer(len(np.arange(f_min, f_max, initial_df / df_inc_factor ** i))).average().save('I')
            Q_st.buffer(len(np.arange(f_min, f_max, initial_df / df_inc_factor ** i))).average().save('Q')
   
    qm = qmm.open_qm(config)
    job = qm.execute(rr_spec)
    job.result_handles.wait_for_all_values()
    I = job.result_handles.get("I").fetch_all()
    Q = job.result_handles.get("Q").fetch_all()

    freq_list = 1e-9*(rr_LO + np.arange(f_min, f_max, initial_df / df_inc_factor ** i))
    zeros = np.where(freq_list == 0)[0]
    if len(zeros) != 0:
        zero_i = zeros[0]
        freq_list = np.delete(freq_list, zero_i)
    
    sig = I + 1j*Q
    if len(zeros) != 0:
        sig = np.delete(sig, zero_i)
    
    e_delay = elec_delay_ns[str(rr_no)]
    p_offset = phase_offset_rad[str(rr_no)]
    if check_e_delay:
        e_delay = elec_delay_ns[str(rr_no)]
        p_offset = phase_offset_rad[str(rr_no)]
        # e_delay = 281
    sig_corrected = sig*np.exp(1j*2*np.pi*freq_list*e_delay + 1j*p_offset)
    phase = np.angle(sig_corrected)
    real = np.real(sig_corrected)
    imag = np.imag(sig_corrected)
    
    abs_sig = np.abs(sig_corrected)
    
    f_min_i = np.argmin(abs_sig)
    f_min = freq_list[f_min_i]
    
    half_max = max(abs_sig) / 2.
    d = abs_sig - half_max
    indexes = np.where(d > 0)[0]
    fwhm = freq_list[indexes[-1]] - freq_list[indexes[0]]
    
    f_min = f_min - fwhm / 2
    f_max = f_min + fwhm
    print("New Range: [{}, {}]".format(f_min, f_max))
    print("Minima: ",f_min)
    
    if(i == 1):
        freq_list_bw = freq_list
    
    plt.plot(freq_list, abs_sig)
    plt.show()
    

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
    return properties
# Calculating the fit parameters and properties
properties = ext_BW(freq_list_bw, real)
print("Total Bandwidth: ", properties['bandwidth'])
print("Internal Bandwidth: ", properties['internal_bandwidth'])
print("External Bandwidth: ", properties['external_bandwidth'])










from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv2 import *
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt


def calculate_FWHM(X, Y):
    half_max = np.min(Y) / 2.0
    left_pts = np.where(Y < half_max)[0]
    right_pts = np.where(Y[::-1] < half_max)[0]
    FWHM = X[-right_pts[0] + 1] - X[left_pts[-1] + 1]
    return FWHM, X[left_pts[-1] + 1], X[-right_pts[0] + 1]


save_data = False
simulate = False
check_e_delay = False

rr_no = 4
rr = f"rr{rr_no}"
out = adc_mapping[rr]
ro_len = ro_len_clk[str(rr_no)]
rep_rate_clk = 2500
rr_LO = config["elements"][rr]["mixInputs"]["lo_frequency"]

f_a = 10e6
f_b = 70e6
df = 0.05e6

if check_e_delay:
    f_a = -250e6
    f_b = 250e6
    df = 1e6

qmm = QuantumMachinesManager()

for i in range(5):  # iteratively update frequency grid
    freq_list = np.arange(f_a, f_b, df)
    zeros = np.where(freq_list == 0)[0]
    if len(zeros):
        zero_i = zeros[0]

    with program() as rr_spec:
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)

        with for_(n, 0, n < 200, n + 1):
            with for_(f, f_a, f < f_b, f + df):

                update_frequency(rr, f)
                wait(rep_rate_clk - ro_len, rr)
                measure("readout", rr, None,
                        demod.full("integW_cos", I, out),
                        demod.full("integW_minus_sin", Q, out))

                save(I, I_st)
                save(Q, Q_st)

        with stream_processing():
            I_st.buffer(len(freq_list)).average().save('I')
            Q_st.buffer(len(freq_list)).average().save('Q')

    if simulate:
        simulation_config = SimulationConfig(
            duration=200000,
            simulation_interface=LoopbackInterface([("con1", 9, "con1", 1), ("con1", 10, "con1", 2)]))
        job = qmm.simulate(config, rr_spec, simulation_config)
        samples = job.get_simulated_samples()
        samples.con1.plot()
        raise Halted()

    qm = qmm.open_qm(config)
    job = qm.execute(rr_spec)
    job.result_handles.wait_for_all_values()
    I = job.result_handles.get("I").fetch_all()
    Q = job.result_handles.get("Q").fetch_all()

    freq_list = 1e-9*(rr_LO + freq_list)
    sig = I + 1j*Q
    if len(zeros):
        freq_list = np.delete(freq_list, zero_i)
        sig = np.delete(sig, zero_i)

    e_delay = elec_delay_ns[str(rr_no)]
    p_offset = phase_offset_rad[str(rr_no)]
    sig_corrected = sig*np.exp(1j*2*np.pi*freq_list*e_delay + 1j*p_offset)
    f_res_i = np.argmin(np.abs(sig))
    f_res = freq_list[f_res_i]

    FWHM, f_a, f_b = calculate_FWHM(freq_list, np.abs(sig))
    df = FWHM / len(freq_list)

print('Minima after 5th iteration:', (f_res, np.min(np.abs(sig))))
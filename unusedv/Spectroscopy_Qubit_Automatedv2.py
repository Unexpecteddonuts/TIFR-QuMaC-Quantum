from qm.qua import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv4 import * # v2 for W lab
from scipy import signal as sgn
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

test = True
if test is True:
    q_no = 2
    
n_avgs = 1e2 # SNR dependent
N_points = 2e3 # No.of points in a sweep window
anharm = 220 # MHz
qmm = QuantumMachinesManager(qop_ip) # Replace with qm_ip for W lab


def q_spectro(q_no, avgs, f_min, f_max, points, q_amp):

    qe = f"q{q_no}"
    rr = f"rr{q_no}"
    out = adc_mapping[rr]
    f_LO = q_LO[f"{q_no}"]

    df = (f_max-f_min)//points
    freqs = np.arange(f_min, f_max, df)

    # freqs = np.round_(np.linspace(f_min, f_max, N_points))
    # df = freqs[1] - freqs[0]

    with program() as qubit_spec:
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)

        with for_(n, 0, n < avgs, n + 1):
            with for_(f, f_min, f < f_max, f + df):
                wait(20000, qe)
                update_frequency(qe, f)
                play("const"*amp(q_amp), qe, duration=20000)
                align(rr,qe)
                measure("readout", rr, None,
                        demod.full("integW_cos", I, out),
                        demod.full("integW_minus_sin", Q, out))
                
                # measure_macro(qe, rr, out, I, Q, pi_12=True)
                save(I, I_st)
                save(Q, Q_st)

        with stream_processing():
            I_st.buffer(len(freqs)).average().save('I')
            Q_st.buffer(len(freqs)).average().save('Q')
    
    return qubit_spec, f_LO, freqs


def run_spectro(params):    

    qua_prog, f_LO, freqs = q_spectro(*params)
    
    qm = qmm.open_qm(config)
    job = qm.execute(qua_prog)
    res_handles = job.result_handles
    I_handle = job.result_handles.get("I")
    Q_handle = job.result_handles.get("Q")
    #job.result_handles.wait_for_all_values()

    plt.figure(1)
    I_handle.wait_for_values(1)
    Q_handle.wait_for_values(1)
    while res_handles.is_processing():

        I = I_handle.fetch_all()
        Q = Q_handle.fetch_all()
        signal = I + 1j * Q
        phase = np.angle(signal)
        amplitude = np.abs(signal)
        sig = amplitude #phase

        mean = np.mean(sig)
        max_sig = np.max(sig)
        min_sig = np.min(sig)

        if (max_sig-mean) < (mean-min_sig):
            sig -= mean
            sig *= -1

        plt.clf()
        plt.plot(1e-6*(freqs), sig, label="Signal")
        plt.xlabel("Frequency offset from LO (MHz)")
        plt.grid()
        plt.legend()
        plt.pause(1)

    I = I_handle.fetch_all()
    Q = Q_handle.fetch_all()
    qm.close()

    signal = I + 1j * Q
    phase = np.angle(signal)
    amplitude = np.abs(signal)
    sig = amplitude #phase

    mean = np.mean(sig)
    max_sig = np.max(sig)
    min_sig = np.min(sig)

    if (max_sig-mean) < (mean-min_sig):
        sig -= mean
        sig *= -1

    plt.clf()
    plt.title("Spectroscopy")
    plt.plot(1e-6*(freqs), sig, label="Signal")
    plt.xlabel("Frequency offset from LO (MHz)")
    plt.grid()
    plt.legend()
    plt.close(1)
    
    freq_array = 1e-9*(f_LO + freqs)
    return sig, freq_array


def broad_sweep(q_no, avgs=n_avgs, f_min_MHz=-350, f_max_MHz=350,
                points=N_points, q_amp=0.4, min_anharm_MHz=anharm-10):
    
    # avgs = 5e2 but for good SNR - 100, bad SNR - 1000
    f_min = f_min_MHz * u.MHz
    f_max = f_max_MHz * u.MHz
    
    df = (f_max-f_min)//points
    min_zero2by2_det = (min_anharm_MHz/2) * u.MHz # zero2by2 detuning from qubit

    spacing = min_zero2by2_det/df

    spectro_params = (q_no, avgs, f_min, f_max, points, q_amp)
    sig, freqs = run_spectro(spectro_params)
    sig -= min(sig)
    sig /= max(sig)

    sig_peaks, peak_props = sgn.find_peaks(sig, prominence=0.4, width=1, rel_height=0.5, distance=spacing)

    # print(f'The signal peaks are found at the indices {sig_peaks}.')
    # print(f'The peak properties are {peak_props}.')

    widths = peak_props['widths']
    q_index = widths.tolist().index(max(widths))
    q_peak = sig_peaks[q_index]
    zero2by2_index = widths.tolist().index(min(widths))
    zero2by2_peak = sig_peaks[zero2by2_index]

    i = 0
    hlines_array = []
    for key, val in peak_props.items():
        # print(f'i is {i}')
        # print(f'key is {key}')
        # print(f'val is {val}')
        
        if i == 3:
            widths = val
        if i==4:
            hlines_array.append(val)
        if i>4:
            hlines_array.append(freqs[0] + val*df*1e-9) # df is in MHz
            
        i += 1
    hlines_list = (hlines_array[0], hlines_array[1], hlines_array[2])

    peaks = np.array([zero2by2_peak, q_peak])
    plt.figure()
    plt.title('Qubit Spectroscopy with 0-2/2 line')
    plt.plot(freqs, sig)
    plt.plot(freqs[sig_peaks], sig[sig_peaks], "x", color='C1')
    plt.plot(freqs[peaks], sig[peaks], "o", color='C2')
    # plt.vlines(x=freqs[sig_peaks], ymin=baseline, ymax=sig[sig_peaks], colors='C2')
    plt.hlines(y=hlines_list[0], xmin=hlines_list[1], xmax=hlines_list[2], colors='C3')
    plt.xlabel('Frequency sweep (GHz)')
    plt.ylabel('Normalized intensity (a.u.)')
    # plt.show()

    q_freq = freqs[q_peak]
    zero2by2_freq = freqs[zero2by2_peak]
    q_fwhm = hlines_list[2][q_index] - hlines_list[1][q_index]  
    zero2by2_fwhm = hlines_list[2][zero2by2_index] - hlines_list[1][zero2by2_index]  

    print('################################\n'
          'Qubit and 0-2/2 Peak Properties\n'
          '================================')
    print(f'The qubit frequency is {q_freq} GHz.')
    print(f'The FWHM of the qubit peak is {1e3 * q_fwhm} MHz.')
    print(f'The zero2by2 frequency is {zero2by2_freq} GHz.')
    print(f'The FWHM of the zero2by2 peak is {1e3 * zero2by2_fwhm} MHz.')
    print('--------------------------------------------------------------')
    zero2by2_det = zero2by2_freq - q_freq
    anharm = 2*(zero2by2_det)
    print(f'The 0-2/2 line is detuned by {1e3 * zero2by2_det} MHz.')
    print(f'The anharmonicity of the qubit is {1e3 * anharm} MHz.')

    rerun_sweep = -zero2by2_det > (50*1e-3 + min_zero2by2_det*1e-9) # Not within 50 MHz of expected 0-2/2 detuning
    print(f"Do I need to rerun broad sweep? A - {rerun_sweep}")
    return [q_amp, rerun_sweep, points], [q_freq, q_fwhm], [zero2by2_freq, zero2by2_fwhm]


def fine_qubit_sweep(q_no, old_q_freq, old_q_fwhm, old_q_amp, 
                     avgs=n_avgs, points=N_points, into_amp=1/3):
    
    f_LO = q_LO[f"{q_no}"]
    # avgs = 5e2 but for good SNR - 100, bad SNR - 1000
    q_offLO = old_q_freq*u.GHz - f_LO
    f_min = q_offLO - 3*old_q_fwhm*u.GHz
    f_max = q_offLO + 3*old_q_fwhm*u.GHz

    new_df = (f_max-f_min)//points # Only needed for printing purposes 
    q_amp = old_q_amp * into_amp
    
    print('#########################################\n'
          'The fine qubit sweep tells us:')
    print(f'The qubit is offset from LO by {q_offLO}.\n'
          f'f_min, f_max, df = {f_min, f_max, new_df}\n'
          f'The control amplitude is {q_amp}.')

    spectro_params = (q_no, avgs, f_min, f_max, points, q_amp)
    sig, freqs = run_spectro(spectro_params)
    sig -= min(sig)
    sig /= max(sig)

    plt.plot(sig)

    sig_peaks, peak_props = sgn.find_peaks(sig, prominence=0.4, width=1, rel_height=0.5)

    widths = peak_props['widths']
    q_index = widths.tolist().index(max(widths))
    q_peak = sig_peaks[q_index]
    
    print('--------------------------------------------------------------')
    print(f'The qubit peak is found at a frequency offset of {q_peak*new_df*1e-6} MHz.')

    for key, val in peak_props.items():
        vars()[f'q_{key[:-1]}'] = val[q_index]

    baseline = 0*sig[q_peak]

    q_freq = freqs[q_peak]
    q_fwhm = max(widths)*new_df*1e-9 # in GHz

    print(f'The qubit frequency is {q_freq} GHz.')
    print(f'The FWHM of the qubit peak is {1e3 * q_fwhm} MHz.')
    
    plt.figure()
    plt.title('Qubit Fine Spectroscopy')
    plt.plot(freqs, sig)
    plt.plot(freqs[q_peak], sig[q_peak], "o", color="C1")
    plt.vlines(x=freqs[q_peak], ymin=baseline, ymax=sig[q_peak], colors='C2')
    plt.hlines(eval(f'q_width_height'), freqs[0] + eval(f'q_left_ip')*new_df*1e-9, freqs[0] + eval(f'q_right_ip')*new_df*1e-9, colors="C3")
    # plt.show()    

    return q_freq, q_fwhm, q_amp, new_df

def fine_zero2by2_sweep(q_no, old_zero2by2_freq, old_zero2by2_fwhm, old_q_amp, old_df, 
               avgs=n_avgs, into_df=1/10, into_amp=15/16):
    
    f_LO = q_LO[f"{q_no}"]
    # avgs = 5e2 but for good SNR - 100, bad SNR - 1000
    zero2by2_offLO = old_zero2by2_freq*u.GHz - f_LO
    f_min = zero2by2_offLO - 3*old_zero2by2_fwhm*u.GHz
    f_max = zero2by2_offLO + 3*old_zero2by2_fwhm*u.GHz
    new_df = old_df * into_df
    q_amp = old_q_amp * into_amp

    print(f'The 0-2/2 line is offset from LO by {zero2by2_offLO}.\n'
          f'f_min, f_max, df = {f_min, f_max, new_df}\n'
          f'The control amplitude is {q_amp}.')

    spectro_params = (q_no, avgs, f_min, f_max, new_df, q_amp)
    sig, freqs = run_spectro(spectro_params)
    sig -= min(sig)
    sig /= max(sig)

    plt.plot(sig)

    sig_peaks, peak_props = sgn.find_peaks(sig, prominence=0.4, width=1, rel_height=0.5)

    widths = peak_props['widths']
    zero2by2_index = widths.tolist().index(max(widths))
    zero2by2_peak = sig_peaks[zero2by2_index]
    
    print(f'The peak is found at the index: {zero2by2_peak}, at a frequency offset of {zero2by2_peak*new_df*1e-6} MHz.')

    for key, val in peak_props.items():
        vars()[f'zero2by2_{key[:-1]}'] = val[zero2by2_index]


    baseline = 0*sig[zero2by2_peak]

    zero2by2_freq = freqs[zero2by2_peak]
    zero2by2_fwhm = max(widths)*new_df*1e-9 # in GHz

    print(f'The 0-2/2 frequency is {zero2by2_freq} GHz.')
    print(f'The FWHM of the 0-2/2 peak is {1e3 * zero2by2_fwhm} MHz.')
    
    plt.figure()
    plt.title('0-2/2 Fine Spectroscopy')
    plt.plot(freqs, sig)
    plt.plot(freqs[zero2by2_peak], sig[zero2by2_peak], "o", color="C1")
    plt.vlines(x=freqs[zero2by2_peak], ymin=baseline, ymax=sig[zero2by2_peak], colors='C2')
    plt.hlines(eval(f'zero2by2_width_height'), freqs[0] + eval(f'zero2by2_left_ip')*new_df*1e-9, freqs[0] + eval(f'zero2by2_right_ip')*new_df*1e-9, colors="C3")
    # plt.show()    

    return zero2by2_freq, zero2by2_fwhm, q_amp, new_df


def auto_spectro(q_no):

    prog_props, q_props, zero2by2_props = broad_sweep(q_no)

    rerun_sweep = prog_props[1]
    print("------------------\n"
         f"Rerun is {rerun_sweep}\n"
          "------------------")
    while rerun_sweep:
        print("--------------------------------------------------------")
        print("Rerunning broad frequency sweep with more power and finer binning.")
        print("--------------------------------------------------------")
        old_amp = prog_props[0]
        old_points = prog_props[2]
        temp_amp = old_amp*5/4
        new_points = old_points*1.1
        prog_props, q_props, zero2by2_props, rerun_sweep  = broad_sweep(q_no, q_amp=temp_amp, points=new_points)


    old_amp = prog_props[0]
    old_points = prog_props[2]
    old_q_freq = q_props[0]
    old_q_fwhm = q_props[1]
    old_zero2by2_freq = zero2by2_props[0]
    old_zero2by2_fwhm = zero2by2_props[1]

        
    q_fwhm = old_q_fwhm * u.GHz
    print(f'Do we need a rerun_sweep of fine sweep? A = {q_fwhm > 50 * u.kHz}')
    while q_fwhm > (50 * u.kHz):
    # for i in range(1,4):
        q_freq, q_fwhm, q_amp, new_df = fine_qubit_sweep(q_no, old_q_freq, old_q_fwhm, old_amp, old_df)
        old_q_freq = q_freq
        old_q_fwhm = q_fwhm
        old_amp = q_amp
        old_df = new_df
        q_fwhm = q_fwhm * u.GHz


    print(f'The final qubit frequency is {q_freq} GHz.')
    print(f'The final FWHM of the qubit peak is {1e-6 * q_fwhm} MHz.')

    zero2by2_fwhm = old_zero2by2_fwhm * u.GHz
    zero2by2_freq = old_zero2by2_freq

    '''while zero2by2_fwhm > 50 * u.kHz:
    # for i in range(1,4):
        zero2by2_freq, zero2by2_fwhm, q_amp, new_df = fine_zero2by2_sweep(q_no, old_zero2by2_freq, old_zero2by2_fwhm, old_amp, old_df)
        old_zero2by2_freq = zero2by2_freq
        old_zero2by2_fwhm = zero2by2_fwhm
        old_amp = q_amp
        old_df = new_df'''

    # q_fwhm = old_q_fwhm * u.GHz
    # while q_fwhm > 1 * u.MHz:
    #     q_freq, q_fwhm = fine_qubit_sweep(q_no, old_q_freq, old_q_fwhm, old_amp, old_df)

    # print(f'The final qubit frequency is {q_freq} GHz.')
    # print(f'The final FWHM of the qubit peak is {1e3 * q_fwhm} MHz.')
    # print('##################################')

    # zero2by2_fwhm = old_zero2by2_fwhm * u.GHz
    # while zero2by2_fwhm > 1 * u.MHz:
    #     zero2by2_freq, zero2by2_fwhm = fine_zero2by2_sweep(q_no, old_zero2by2_freq, old_zero2by2_fwhm, old_amp, old_df)

    print(f'The final zero2by2 frequency is {zero2by2_freq} GHz.')
    print(f'The final FWHM of the zero2by2 peak is {1e3 * zero2by2_fwhm} MHz.')
    print('##################################')

    zero2by2_det = zero2by2_freq - q_freq
    anharm = 2*zero2by2_det

    print(f'The 0-2/2 line is detuned by {1e3 * zero2by2_det} MHz.')
    print(f'The final anharmonicity of the qubit is {1e3 * anharm} MHz.')

    return [q_freq, q_fwhm], [zero2by2_freq, zero2by2_fwhm], anharm

if test is True:
    auto_spectro(q_no)
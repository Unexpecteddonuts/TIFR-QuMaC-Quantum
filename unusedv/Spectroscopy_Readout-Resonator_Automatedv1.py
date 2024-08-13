import json
from qm import SimulationConfig
from qm.qua import *
from qm import LoopbackInterface
from qm.QuantumMachinesManager import QuantumMachinesManager
from configuration_4qubitsv2 import *
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
save_data = False
simulate =True

###################
# The QUA program #
###################
res_freq = []

f_min = -10e6
f_max = 70e6
df = 0.1e6

def rr_spectro(res_freq, rr_no, f_min, f_max, df, check_e_delay=False ):
    
    out = adc_mapping[rr]
    ro_len = ro_len_clk[str(rr_no)]
    rep_rate_clk = 2500
    rr_LO = config["elements"][rr]["mixInputs"]["lo_frequency"]
    rr = f"rr{rr_no}"
    
    if check_e_delay:
        f_min = -250e6
        f_max = 250e6
        df = 1e6

    freq_list = np.arange(f_min, f_max, df)
    zero_i = np.where(freq_list == 0)[0][0]

    with program() as rr_spec:
        n = declare(int)
        I = declare(fixed)
        I_st = declare_stream()
        Q = declare(fixed)
        Q_st = declare_stream()
        f = declare(int)

        with for_(n, 0, n < 100, n + 1):
            with for_(f, f_min, f < f_max, f + df):

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

    ######################################
    # Open Communication with the Server #
    ######################################
    qmm = QuantumMachinesManager(opx_ip)
    
    #############
    # execution #
    #############

    qm = qmm.open_qm(config)
    job = qm.execute(rr_spec)
    job.result_handles.wait_for_all_values()
    I = job.result_handles.get("I").fetch_all()
    Q = job.result_handles.get("Q").fetch_all()

    ############
    # analysis #
    ############
    freq_list = 1e-9*(rr_LO + freq_list)
    sig = I + 1j*Q
    freq_list = np.delete(freq_list, zero_i)
    sig = np.delete(sig, zero_i)

    e_delay = elec_delay_ns[str(rr_no)]
    p_offset =phase_offset_rad[str(rr_no)]
    sig_corrected = sig*np.exp(1j*2*np.pi*freq_list*e_delay + 1j*p_offset)
    phase = np.angle(sig_corrected)
    real = np.real(sig_corrected)
    f_res_i = np.argmin(abs(phase))
    f_res = freq_list[f_res_i]
    res_freq['rr_'+str(rr_no)] = f_res

    data_p = np.column_stack((freq_list, phase, real, np.abs(sig)))
    np.savetxt(f'rr_{rr_no}_spec', data_p, delimiter=",", header="Frequency (GHz), Phase (rad), Real, Magnitude")
    
    return res_freq


for rr_no in range(1,4):
    rr_spectro(res_freq, rr_no, f_min, f_max, df, check_e_delay=False )
    with open('rrfreq.json', 'w') as jsonfile:
     json.dump(res_freq, jsonfile)
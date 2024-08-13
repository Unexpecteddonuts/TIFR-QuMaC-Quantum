from config_builder import *
from helper_functionsv2 import *
from qualang_tools.units import unit
from instrumentlib import RhodeandSchwarz_SGS100A
import matplotlib
import json
import values
values.modify_config_values()

u = unit()
matplotlib.use('Qt5Agg')

opx_ip = "192.168.0.106"
Line4_Line3 = False
TOF_testing = False
initialize_RF_sources = True
path = r"D:\Experiments\2023-05-20 Ringv1_4Q_INDIQ05_Run2"
ExpName = "2023-05-20 Ringv1_4Q_INDIQ05_Run2"

# ====================== DAC Mapping ========================================
dac_mapping = {"q1": [1, [1, 2]], "q2": [1, [3, 4]], "q3": [2, [1, 2]], "q4": [2, [3, 4]],
               "rr1": [1, [7, 8]], "rr2": [1, [9, 10]], "rr3": [2, [7, 8]], "rr4": [2, [9, 10]],
               "cr_c2t1": [1, [3, 4]], "cr_c4t1": [2, [3, 4]], "cr_c2t3": [1, [3, 4]], "cr_c4t3": [2, [3, 4]]}

adc_mapping = {"rr1": "out1", "rr2": "out2", "rr3": "out1", "rr4": "out2"}

# =======================Frequencies==========================================

qe_list = ["q1", "q2", "q3", "q4",
           "rr1", "rr2", "rr3", "rr4",
           "cr_c2t1", "cr_c4t1", "cr_c2t3", "cr_c4t3",
           "q12_1", "q12_2", "q12_3", "q12_4"
           ]

CrossKerr = {"1": 0.115 * u.MHz,              # Calculate and update later using simulation/ conditional Ramsey
             "2": 0.0975 * u.MHz,
             "3": 0.1095 * u.MHz,
             "4": 0.1255 * u.MHz,
             }

q_LO = {"1": 4.9 * u.GHz,
        "2": 4.9 * u.GHz,
        "3": 4.9 * u.GHz,
        "4": 4.9 * u.GHz,
        }

q_IF = {"1": 60.98 * u.MHz + CrossKerr["1"],    #Ring 1
        "2": 88.6733 * u.MHz + CrossKerr["2"],   #Ring 11
        "3": -46.83847 * u.MHz + CrossKerr["3"],   #Ring 5
        "4": 104.23544 * u.MHz + CrossKerr["4"],  #Ring 7
        }

rr_LO = {"1": 7.414055 * u.GHz,
         "2": 7.43615 * u.GHz, #7.43002
         "3": 7.362655 * u.GHz,
         "4": 7.362655 * u.GHz,
         }

rr_IF = {"1": 20 * u.MHz,
         "2": 20 * u.MHz,
         "3": 40 * u.MHz, #39.49
         "4": 36.425 * u.MHz, # 36.985  * u.MHz
         }

# ======================Readout Parameters===========================================

tof = {"1": 300,
       "2": 300,
       "3": 300,
       "4": 300,
       }

ro_len_clk = {"1": 500,
              "2": 500,
              "3": 500,
              "4": 2*500,
              }

ro_amp = {"1": 0.216,
          "2": 0.3,
          "3": 0.192, # 0.4,
          "4": 0.04,           # earlier 0.03, got slightly better with 1.5 times.
          }

# For TOF Testing
if TOF_testing:

    ro_amp = {"1": 1.0,
              "2": 1.0,
              "3": 1.0,
              "4": 1.0,
              }

    rr_IF = {"1": 12 * u.MHz,
             "2": 12 * u.MHz,
             "3": 12 * u.MHz,
             "4": 12 * u.MHz,
             }

integ_len_clk = {"1": 500,
                 "2": 500,
                 "3": 500,
                 "4": 2*500,
                 }

optimal_readout_phase = {"rr1": 123.4*(np.pi/180),
                         "rr2": -188*(-np.pi/180),
                         "rr3": -7.34*(np.pi/180),
                         "rr4": (33.6)*(np.pi/180), #33.6 - 74.7
                         }

demarcations = {"1": 5.475e-05,
                "2": -1.997e-04,
                "3": 5.109e-05,
                "4": -1.190e-05,
                }

elec_delay_ns = {"1": 273.45,
                "2": 280,
                "3": 277.5,
                "4": 277.5,
                }

phase_offset_rad = {"1": -3.269,
                "2": 4.494,
                "3": 2.14,
                "4": 3.19 - 1.6,
                }

# ==================Control Parameters===============================================

pi_rise_grft_ns = 10
pi_len_ns = {"1": 52,
             "2": 52,
             "3": 72,
             "4": 52,
             }

piby2_rise_grft_ns = 10
piby2_len_ns = {"1": 52,
                "2": 52,
                "3": 72,
                "4": 52,
                }

cr_tail_ns = 16
cr_len_ns = {"cr_c2t1": 336,                 # CR gate length except 16ns rise and 16ns fall
             "cr_c4t1": 412,
             "cr_c2t3": 188,
             "cr_c4t3": 576,
             }

# cr_amp = {"cr_c2t1": -0.3, "cr_ac_c2t1": 0.03,
#           "cr_c4t1": -0.05, "cr_ac_c4t1": 0.03, # 0.25/2  #0.25/3 tested, good RB and echo
#           "cr_c2t3": 0.6, "cr_ac_c2t3": 0.03,
#           "cr_c4t3": -0.08, "cr_ac_c4t3": 0.03}  # RB Yes

cr_amp = {"cr_c2t1": -0.3, "cr_ac_c2t1": 0.03,
          "cr_c4t1": 0.05, "cr_ac_c4t1": 0.03,
          "cr_c2t3": -0.6, "cr_ac_c2t3": 0.03,
          "cr_c4t3": 0.08, "cr_ac_c4t3": 0.03}  #CN Yes

cr_phase = {"cr_c2t1":   0.4474474474474474, "cr_ac_c2t1": 0.13713713713713713,
            "cr_c4t1":   0.16516516516516516, "cr_ac_c4t1": -0.32132132132132135,
            "cr_c2t3":   0.3883883883883884, "cr_ac_c2t3": 0.08408408408408413,
            "cr_c4t3":   0.4734734734734734, "cr_ac_c4t3": 0.03003003003003002}

calib_vals = {"1": {"amin": 0.46, "amax": 0.52, "da": 0.0005, "n_pulses": 13},
              "2": {"amin": 0.295, "amax": 0.318, "da": 0.0001, "n_pulses": 17},
              "3": {"amin": 0.35, "amax": 0.50, "da": 0.001, "n_pulses": 7},
              "4": {"amin": 0.12, "amax": 0.126, "da": 0.00005, "n_pulses": 17}
              }


#ojnoirnfopwerifnoergnotingitngoeirgoeingingoringoinoihnotinhoihntoihnttihntohithn


amp_scale = {"1": {"X180": 0.48983491745, "Y180":  0.4896248124, "X90": 0.243341670, "Y90": 0.2433266633},
             "2": {"X180": 0.30521710, "Y180": 0.30507253626, "X90":  0.15164707353, "Y90":  0.15164707353},
             "3": {"X180": 0.8499399, "Y180": 0.8499399, "X90": 0.42496995, "Y90": 0.42496995},  #0.5981540770
             "4": {"X180": 0.123298649, "Y180": 0.123298649, "X90": 0.06149624812, "Y90": 0.06149624812}, #0.1233579289
             }

# =================== Dictionaries for 1-2 transition =================================

q12_IF = {"1": (60.98 - 303.8) * u.MHz,    #Ring 1
          "2": (88.6733 - 309.8) * u.MHz,   #Ring 11
          "3": (-46.83847 - 311.2) * u.MHz,   #Ring 5
          "4": -204.595 * u.MHz,  #Ring 7
         }

piby2_12_len_ns = {"1": 52,
                "2": 52,
                "3": 72,
                "4": 52,
                }

pi_12_len_ns = {"1": 52,
                "2": 52,
                "3": 72,
                "4": 52,
                }

amp_12_scale = {"1": {"X180": 0.48983491745, "Y180":  0.4896248124, "X90": 0.243341670, "Y90": 0.2433266633},
                "2": {"X180": 0.30521710, "Y180": 0.30507253626, "X90":  0.15164707353, "Y90":  0.15164707353},
                "3": {"X180": 0.86133, "Y180": 0.45565282, "X90": 0.227151075, "Y90": 0.227151075},  #0.5981540770
                "4": {"X180": 0.1104652326, "Y180": 0.1104652326, "X90": 0.05494747, "Y90": 0.05494747}, #0.1233579289
                }

# =================== Debugging Mode =================================================
if Line4_Line3:

    CrossKerr["4"], CrossKerr["3"] = CrossKerr["3"], CrossKerr["4"]

    q_LO["4"], q_LO["3"] = q_LO["3"], q_LO["4"]
    rr_LO["4"], rr_LO["3"] = rr_LO["3"], rr_LO["4"]

    q_IF["4"], q_IF["3"] = q_IF["3"], q_IF["4"]
    rr_IF["4"], rr_IF["3"] = rr_IF["3"], rr_IF["4"]

    ro_amp["4"], ro_amp["3"] = ro_amp["3"], ro_amp["4"]
    amp_scale["4"], amp_scale["3"] = amp_scale["3"], amp_scale["4"]
    tof["4"], tof["3"] = tof["3"], tof["4"]

    optimal_readout_phase["rr4"], optimal_readout_phase["rr3"] = optimal_readout_phase["rr3"], optimal_readout_phase["rr4"]

# ==================Mixer Parameters ===============================================
# =================== DC Offsets ===============================================
with open('/home/kalibot/QASMtoQUA/QASMtoExec/fourqubitv2/Calibrations/dc_offsets.json','r') as f:
    dc_offsets = json.load(f)

# ================ List of Mixers ===============================================
mixers = {}
for qe in qe_list:
    mixers[qe] = "mixer_" + qe

# ================ IQ Imbalance Correction Matrices ===============================================
with open('/home/kalibot/QASMtoQUA/QASMtoExec/fourqubitv2/Calibrations/iq_imbalance.json','r') as f:
    iq_imbalance = json.load(f)

mixer_corrections = {}
for qe in qe_list:
    a = iq_imbalance[qe]["a"]
    p = iq_imbalance[qe]["p"]
    mixer_corrections[qe] = IQ_imbalance(a, p)

# ============== ADC Offets ===============================================
with open('/home/kalibot/QASMtoQUA/QASMtoExec/fourqubitv2/Calibrations/adc_offsets.json','r') as f:
    adc_offsets = json.load(f)

# ======================= Setting up RF sources with above parameters =========================

def setup_RF_Source(rf_ip, rf_f, rf_p):

    rf = RhodeandSchwarz_SGS100A(rf_ip)
    rf.set_freq_Hz(rf_f)
    rf.set_power_dB(rf_p)
    rf.output_on()
    rf.close_connection()


rf_rr1_ip = "TCPIP0::192.168.0.121::inst0::INSTR"
rf_rr2_ip = "TCPIP0::192.168.0.108::inst0::INSTR"
rf_rr34_ip = "TCPIP0::192.168.0.200::inst0::INSTR"

rf_q12_ip = "TCPIP0::192.168.0.104::inst0::INSTR"
rf_q34_ip = "TCPIP0::192.168.0.103::inst0::INSTR"

rr_LO_dBm = 19
q_LO_dBm = 17

if initialize_RF_sources:
    setup_RF_Source(rf_rr1_ip, rr_LO["1"], rr_LO_dBm)
    setup_RF_Source(rf_rr2_ip, rr_LO["2"], rr_LO_dBm)
    setup_RF_Source(rf_rr34_ip, rr_LO["3"], rr_LO_dBm+3)
    setup_RF_Source(rf_q12_ip, q_LO["1"], q_LO_dBm)

# ======================== QM Config Starts Here =========================================

config = {

    "version": 1,
}

config = config_add_controller(config, 1, dc_offsets, adc_offsets)
config = config_add_controller(config, 2, dc_offsets, adc_offsets)
config = config_add_controller(config, 3, dc_offsets, adc_offsets)

config = config_add_common_elements(config)

n_qubits = 4
for i in range(1, n_qubits+1):
    q_no = i
    rr_no = i

    config = config_add_elements_q_rr(config, q_no, rr_no, dac_mapping, q_LO, q_IF, rr_LO, rr_IF, pi_len_ns, piby2_len_ns, pi_rise_grft_ns, amp_scale,
                                      mixers, mixer_corrections, ro_amp, ro_len_clk, tof, integ_len_clk, optimal_readout_phase, smearing=0)

config = config_add_rise_fall(config, cr_tail_ns)

# Check CR config elements
config = config_add_crgate(config, 2, 1, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)
config = config_add_crgate(config, 4, 1, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)
config = config_add_crgate(config, 2, 3, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)
config = config_add_crgate(config, 4, 3, dac_mapping, q_LO, q_IF, mixers, mixer_corrections)


#Add quantum elements for 1-2 transition
for i in range(1, n_qubits+1):

    q_no = i
    config = config_add_q12(config, q_no, dac_mapping, q_LO, q12_IF,
                   pi_12_len_ns, piby2_12_len_ns, pi_rise_grft_ns, amp_12_scale, mixers, mixer_corrections)
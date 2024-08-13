from helper_functionsv2 import *

def config_add_controller(config, con_no, dc_offsets, adc_offsets):

    if not "controllers" in config: config["controllers"] = {}

    config["controllers"][f"con{con_no}"] = {}
    config["controllers"][f"con{con_no}"]["type"] = "opx1"
    config["controllers"][f"con{con_no}"]["analog_outputs"] = {}

    for i in range(1, 11):
        config["controllers"][f"con{con_no}"]["analog_outputs"][i] = {"offset" : dc_offsets[f"con{con_no}"][f"{i}"]}

    config["controllers"][f"con{con_no}"]["digital_outputs"] = {}
    config["controllers"][f"con{con_no}"]["analog_inputs"] = {}

    for i in range(1, 3):
        config["controllers"][f"con{con_no}"]["analog_inputs"][i] = {"offset": adc_offsets[f"con{con_no}"][f"{i}"], 'gain_db': 20}

    return config

def config_add_common_elements(config):
    config["pulses"] = {"zero": {"operation": "control",
                                 "length": 100,
                                 "waveforms": {"I": "zero_wf", "Q": "zero_wf"},
                                 },

                        "const_pulse": {
                            "operation": "control",
                            "length": 100,
                            "waveforms": {"I": "const_wf", "Q": "zero_wf"},
                        },
                        }

    config["waveforms"] = {"zero_wf": {"type": "constant", "sample": 0.0},
                           "const_wf": {"type": "constant", "sample": 0.4},
                           }

    config["digital_waveforms"] = {"ON": {"samples": [(1, 0)]}}

    return config

def config_add_elements_q_rr(config, q_no, rr_no, dac_mapping, q_LO, q_IF, rr_LO, rr_IF, pi_len_ns, piby2_len_ns, pi_rise_grft_ns, amp_scale,
                     mixers, mixer_corrections, ro_amp, ro_len_clk, tof, integ_len_clk, optimal_readout_phase, smearing=0):
    '''
    Given parameters, this function adds the qubit, readout resonator, and all the standard pulses and waveforms to the config.
    '''

    ################ ADD ELEMENTS ################
    # Add qubit
    if not "elements" in config: config["elements"] = {}
    config['elements'][f"q{q_no}"] = {
        "mixInputs": {
            "I": (f"con{dac_mapping[f'q{q_no}'][0]}", dac_mapping[f'q{q_no}'][1][0]),
            "Q": (f"con{dac_mapping[f'q{q_no}'][0]}", dac_mapping[f'q{q_no}'][1][1]),
            "lo_frequency": q_LO[f"{q_no}"],
            "mixer": mixers[f"q{q_no}"],
        },
        "intermediate_frequency": q_IF[f"{q_no}"],
        "operations": {
            "const": "const_pulse",
            "grft": f"q{q_no}_grft",
            "X180": f"q{q_no}_X180",
            "X90": f"q{q_no}_X90",
            "mX90": f"q{q_no}_mX90",
            "Y180": f"q{q_no}_Y180",
            "Y90": f"q{q_no}_Y90",
            "mY90": f"q{q_no}_mY90",
            "I": "zero",
        },
    }

    # add corresponding rr
    config['elements'][f"rr{rr_no}"] = {
        "mixInputs": {
            "I": (f"con{dac_mapping[f'rr{rr_no}'][0]}", dac_mapping[f'rr{rr_no}'][1][0]),
            "Q": (f"con{dac_mapping[f'rr{rr_no}'][0]}", dac_mapping[f'rr{rr_no}'][1][1]),
            "lo_frequency": rr_LO[f"{rr_no}"],
            "mixer": mixers[f"rr{rr_no}"],
        },
        "intermediate_frequency": rr_IF[f"{rr_no}"],
        "outputs": {
            "out1": (f"con{dac_mapping[f'rr{rr_no}'][0]}", 1),
            "out2": (f"con{dac_mapping[f'rr{rr_no}'][0]}", 2),
        },
        "time_of_flight": tof[f"{rr_no}"],
        "smearing": smearing,
        "operations": {
            "const": f"const_pulse",
            "readout": f"q{rr_no}_ro_pulse",
        }
    }

    ####### ADD CORRESPONDING PULSES #########
    if not "pulses" in config: config["pulses"] = {}
    if not "waveforms" in config: config["waveforms"] = {}
    # add all control pulses
    for key in config["elements"][f"q{q_no}"]["operations"]:
        pul = config["elements"][f"q{q_no}"]["operations"][key]
        # print(pul)
        if pul == "const_pulse" or pul == "zero":
            continue

        pul_name = pul.split('_')[1]
        if "m" in pul: pul_name = pul.split('_')[1][1:]


        config["pulses"][pul] = {
            "operation": "control",
            "length": pi_len_ns[f"{q_no}"],
            "waveforms": {"I": "zero_wf", "Q": "zero_wf"},
        }

        if 'grft' in pul:
            config["pulses"][pul]["waveforms"]["I"] = f"{pul}_I_wf"
            config["waveforms"][f"{pul}_I_wf"] = {"type": "arbitrary",
                                                  "samples": grft_arr_gen((pi_len_ns[f"{q_no}"], pi_rise_grft_ns), 1)}

        if "X" in pul:
            config["pulses"][pul]["waveforms"]["I"] = f"{pul}_I_wf"
            # add waveform
            if "m" in pul:
                config["waveforms"][f"{pul}_I_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [-amp_scale[f"{q_no}"][pul_name]])}
            else:
                config["waveforms"][f"{pul}_I_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [amp_scale[f"{q_no}"][pul_name]])}

        if "Y" in pul:
            config["pulses"][pul]["waveforms"]["Q"] = f"{pul}_Q_wf"
            # add waveform
            if "m" in pul:
                config["waveforms"][f"{pul}_Q_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [-amp_scale[f"{q_no}"][pul_name]])}
            else:
                config["waveforms"][f"{pul}_Q_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [amp_scale[f"{q_no}"][pul_name]])}

        if "90" in pul:
            config["pulses"][pul]["length"] = piby2_len_ns[f"{q_no}"]

    # add readout pulse
    config["pulses"][f"q{q_no}_ro_pulse"] = {
        "operation": "measurement",
        "length": ro_len_clk[f"{q_no}"] * 4,
        "waveforms": {"I": f"q{q_no}_ro_wf", "Q": "zero_wf"},
        "integration_weights": {
            "integW_cos": f"integW_cos_rr{q_no}",
            "integW_sin": f"integW_sin_rr{q_no}",
            "integW_minus_sin": f"integW_minus_sin_rr{q_no}"
        },
        "digital_marker": "ON",
    }

    ro_pulse_square = False
    if ro_pulse_square:
        config["waveforms"][f"q{rr_no}_ro_wf"] = {"type": "constant",
                                                  "sample": 0.4 * ro_amp[str(rr_no)]}

    else:
        config["waveforms"][f"q{rr_no}_ro_wf"] = {"type": "arbitrary",
                                                  "samples": grft_arr_gen((ro_len_clk[str(rr_no)] * 4, 200),
                                                  [ro_amp[str(rr_no)]])}

    # add integration weights
    if not "integration_weights" in config: config["integration_weights"] = {}

    config["integration_weights"][f"integW_cos_rr{rr_no}"] = {}
    config["integration_weights"][f"integW_sin_rr{rr_no}"] = {}
    config["integration_weights"][f"integW_minus_sin_rr{rr_no}"] = {}


    config["integration_weights"][f"integW_cos_rr{rr_no}"]["cosine"] = [
        (np.cos(optimal_readout_phase[f"rr{rr_no}"]), integ_len_clk[f"{rr_no}"] * 4)]
    config["integration_weights"][f"integW_cos_rr{rr_no}"]["sine"] = [
        (-np.sin(optimal_readout_phase[f"rr{rr_no}"]), integ_len_clk[f"{rr_no}"] * 4)]

    config["integration_weights"][f"integW_sin_rr{rr_no}"]["cosine"] = [
        (np.sin(optimal_readout_phase[f"rr{rr_no}"]), integ_len_clk[f"{rr_no}"] * 4)]
    config["integration_weights"][f"integW_sin_rr{rr_no}"]["sine"] = [
        (np.cos(optimal_readout_phase[f"rr{rr_no}"]), integ_len_clk[f"{rr_no}"] * 4)]

    config["integration_weights"][f"integW_minus_sin_rr{rr_no}"]["cosine"] = [
        (-np.sin(optimal_readout_phase[f"rr{rr_no}"]), integ_len_clk[f"{rr_no}"] * 4)]
    config["integration_weights"][f"integW_minus_sin_rr{rr_no}"]["sine"] = [
        (-np.cos(optimal_readout_phase[f"rr{rr_no}"]), integ_len_clk[f"{rr_no}"] * 4)]

    # add mixers
    if not "mixers" in config: config["mixers"] = {}
    config["mixers"][f"mixer_q{q_no}"] = [{"intermediate_frequency": q_IF[f"{q_no}"], "lo_frequency": q_LO[f"{q_no}"],
                                           "correction": mixer_corrections[f"q{q_no}"]}]
    config["mixers"][f"mixer_rr{q_no}"] = [
        {"intermediate_frequency": rr_IF[f"{rr_no}"], "lo_frequency": rr_LO[f"{rr_no}"],
         "correction": mixer_corrections[f"rr{rr_no}"]}]

    return config

def config_add_rise_fall(config, cr_tail_ns):

    config["pulses"]["rise_pulse"] = {
            "operation": "control",
            "length": cr_tail_ns,
            "waveforms": {"I": "rise_wf", "Q": "zero_wf"},
        }

    config["pulses"]["fall_pulse"] = {
        "operation": "control",
        "length": cr_tail_ns,
        "waveforms": {"I": "fall_wf", "Q": "zero_wf"},
    }

    config["waveforms"]["rise_wf"] = {"type": "arbitrary", "samples": rise_arr(cr_tail_ns)}
    config["waveforms"]["fall_wf"] = {"type": "arbitrary", "samples": fall_arr(cr_tail_ns)}

    return config


def config_add_crgate(config, control, target, dac_mapping, q_LO, q_IF, mixers, mixer_corrections):

    qe = f"cr_c{control}t{target}"
    qe_ac = f"cr_ac_c{control}t{target}"

    config['elements'][qe] = {
        "mixInputs": {
            "I": (f"con{dac_mapping[f'q{control}'][0]}", dac_mapping[qe][1][0]),
            "Q": (f"con{dac_mapping[f'q{control}'][0]}", dac_mapping[qe][1][1]),
            "lo_frequency": q_LO[f"{control}"],
            "mixer": mixers[qe],
        },
        "intermediate_frequency": q_IF[f"{target}"],
        "operations": {
            "const": "const_pulse",
            "rise": "rise_pulse",
            "fall": "fall_pulse",
        },
    }

    config['elements'][qe_ac] = {
        "mixInputs": {
            "I": (f"con{dac_mapping[f'q{target}'][0]}", dac_mapping[f'q{target}'][1][0]),
            "Q": (f"con{dac_mapping[f'q{target}'][0]}", dac_mapping[f'q{target}'][1][0]),
            "lo_frequency": q_LO[f"{control}"],
            "mixer": mixers[f'q{target}'],
        },
        "intermediate_frequency": q_IF[f"{target}"],
        "operations": {
            "const": "const_pulse",
            "rise": "rise_pulse",
            "fall": "fall_pulse",
        },
    }

    config["mixers"][mixers[qe]] = [{"intermediate_frequency": q_IF[f"{target}"], "lo_frequency": q_LO[f"{control}"], "correction": mixer_corrections[qe]}]

    return config


def config_add_q12(config, q_no, dac_mapping, q_LO, q12_IF,
                   pi_12_len_ns, piby2_12_len_ns, pi_rise_grft_ns, amp_12_scale, mixers, mixer_corrections):

    config['elements'][f"q12_{q_no}"] = {
        "mixInputs": {
            "I": (f"con{dac_mapping[f'q{q_no}'][0]}", dac_mapping[f'q{q_no}'][1][0]),
            "Q": (f"con{dac_mapping[f'q{q_no}'][0]}", dac_mapping[f'q{q_no}'][1][1]),
            "lo_frequency": q_LO[f"{q_no}"],
            "mixer": mixers[f"q12_{q_no}"],
        },
        "intermediate_frequency": q12_IF[f"{q_no}"],
        "operations": {
            "const": "const_pulse",
            "grft": f"q12_{q_no}_grft",
            "X180": f"q12_{q_no}_X180",
            "X90": f"q12_{q_no}_X90",
            "mX90": f"q12_{q_no}_mX90",
            "Y180": f"q12_{q_no}_Y180",
            "Y90": f"q12_{q_no}_Y90",
            "mY90": f"q12_{q_no}_mY90",
            "I": "zero",
        },
    }

    for key in config["elements"][f"q12_{q_no}"]["operations"]:
        pul = config["elements"][f"q12_{q_no}"]["operations"][key]
        # print(pul)
        if pul == "const_pulse" or pul == "zero":
            continue

        pul_name = pul.split('_')[2]
        if "m" in pul: pul_name = pul.split('_')[2][1:]


        config["pulses"][pul] = {
            "operation": "control",
            "length": pi_12_len_ns[f"{q_no}"],
            "waveforms": {"I": "zero_wf", "Q": "zero_wf"},
        }

        if 'grft' in pul:
            config["pulses"][pul]["waveforms"]["I"] = f"{pul}_I_wf"
            config["waveforms"][f"{pul}_I_wf"] = {"type": "arbitrary",
                                                  "samples": grft_arr_gen((pi_12_len_ns[f"{q_no}"], pi_rise_grft_ns), 1)}

        if "X" in pul:
            config["pulses"][pul]["waveforms"]["I"] = f"{pul}_I_wf"
            # add waveform
            if "m" in pul:
                config["waveforms"][f"{pul}_I_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_12_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [-amp_12_scale[f"{q_no}"][pul_name]])}
            else:
                config["waveforms"][f"{pul}_I_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_12_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [amp_12_scale[f"{q_no}"][pul_name]])}

        if "Y" in pul:
            config["pulses"][pul]["waveforms"]["Q"] = f"{pul}_Q_wf"
            # add waveform
            if "m" in pul:
                config["waveforms"][f"{pul}_Q_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_12_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [-amp_12_scale[f"{q_no}"][pul_name]])}
            else:
                config["waveforms"][f"{pul}_Q_wf"] = {"type": "arbitrary",
                                                      "samples": grft_arr_gen((pi_12_len_ns[f"{q_no}"], pi_rise_grft_ns),
                                                                [amp_12_scale[f"{q_no}"][pul_name]])}

        if "90" in pul:
            config["pulses"][pul]["length"] = piby2_12_len_ns[f"{q_no}"]

        config["mixers"][f"mixer_q12_{q_no}"] = [
            {"intermediate_frequency": q12_IF[f"{q_no}"], "lo_frequency": q_LO[f"{q_no}"],
             "correction": mixer_corrections[f"q12_{q_no}"]}]

    return config
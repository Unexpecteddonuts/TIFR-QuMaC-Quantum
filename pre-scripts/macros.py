from qm.qua import *
from configuration_4qubitsv2 import cr_amp, cr_len_ns


def play_flat_top(qe, a, t = 5):

    play("rise"*amp(a), qe)
    play("const"*amp(a), qe, duration=t)
    play("fall"*amp(a), qe)


def Hadamard(qe):

    # frame_rotation_2pi(-0.5, qe)
    play("Y90", qe)
    play("X180", qe)
    play("Y180", qe)


def ZXby4(qe_cr, qe_c, a=1.0, t=28):

    tby2 = t // 2
    play_flat_top(qe_cr, a, tby2)
    align(qe_c, qe_cr)


def ZXby2_echo_noAC(qe_cr, qe_c, qe_t):

    a = cr_amp[qe_cr]
    t = cr_len_ns[qe_cr]//4

    tby2 = t // 2

    align(qe_cr, qe_c, qe_t)
    play_flat_top(qe_cr, a, tby2)
    align(qe_c, qe_cr)
    wait(4, qe_c)
    play("X180", qe_c)
    align(qe_c, qe_cr)
    play_flat_top(qe_cr, -a, tby2)
    align(qe_c, qe_cr)
    wait(4, qe_c)
    play("X180", qe_c)


def CNOT_macro(qe_c, qe_t):

    c, t = qe_c[-1], qe_t[-1]

    if qe_c == "q2" or qe_c == "q4":

        qe_cr = f"cr_c{c}t{t}"
        ZXby2_echo_noAC(qe_cr, qe_c, qe_t)
        align(qe_cr, qe_t)
        wait(4, qe_t)
        play("X90", qe_t)  #fixed for originally mX90
        wait(4, qe_t)
        align(qe_c, qe_t)
        play("X90", qe_c)
        play("Y90", qe_c)
        play("mX90", qe_c)

    if qe_c == "q1" or qe_c == "q3":

        qe_cr = f"cr_c{t}t{c}"
        align(qe_cr, qe_c, qe_t)
        Hadamard(qe_c)
        Hadamard(qe_t)
        wait(4, qe_c)
        align(qe_cr, qe_c, qe_t)
        CNOT_macro(qe_t, qe_c)
        wait(4, qe_c)
        align(qe_cr, qe_c, qe_t)
        Hadamard(qe_c)
        Hadamard(qe_t)

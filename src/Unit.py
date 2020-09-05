import numpy as np
import scipy.constants as cnt

SI = {
    'um': 10**(-6), 'mm': 0.001, 'cm': 0.01, 'm': 1.0, 'km': 1000.,
    'GHz': 10**9, 'MHz': 10**6, 'kHz': 10**3, 'Hz': 1.0,
    'ns': 10**(-9), 'us': 10**(-6), 'ms': 10**(-3), 's': 1.0,
    'min': 60.0, 'h': 60.0 * 60.0, 'day': 60.0 * 60.0 * 24.0,
    'month': 60.0 * 60.0 * 24.0 * 30,
    'year': 60.0 * 60.0 * 24.0 * 365
}


def convert_SI(val_unit, unit_in="cm", unit_out="mm"):
    val = float(val_unit.split(unit_in)[0])
    return val * SI[unit_in] / SI[unit_out]


def convert(unit_in="cm", unit_out="mm"):
    return SI[unit_in] / SI[unit_out]


def convert_freq_to_time(num=1.0, unit_in="Hz", unit_out="s"):
    freq = num * SI[unit_in] / SI["Hz"]
    time_s = 1 / freq
    return time_s * SI["s"] / SI[unit_out]


def convert_time_to_freq(num=1.0, unit_in="s", unit_out="Hz"):
    time_s = num * SI[unit_in] / SI["s"]
    freq = 1 / time_s
    return freq * SI["Hz"] / SI[unit_out]


def knum_from_freq(num=100.0, unit_in="GHz"):
    freq = num * convert(unit_in, "Hz")
    wave = cnt.c / freq
    knum = 2 * np.pi / wave
    return knum, wave


def knum_from_wave(num=1.0, unit_in="mm"):
    wave = num * convert(unit_in, "m")
    freq = cnt.c / wave
    knum = 2 * np.pi / wave
    return knum, freq


if __name__ == "__main__":
    # freq -> time
    print(convert_freq_to_time(1.0, "Hz", "s"), "s")
    print(convert_freq_to_time(2.0, "Hz", "s"), "s")
    print(convert_freq_to_time(2.0, "kHz", "ms"), "ms")
    print(convert_freq_to_time(10.0, "kHz", "us"), "us")
    # time -> freq
    print(convert_time_to_freq(1.0, "s", "Hz"), "Hz")
    print(convert_time_to_freq(0.5, "s", "Hz"), "Hz")
    print(convert_time_to_freq(100, "ms", "Hz"), "Hz")
    print(convert_time_to_freq(100, "us", "kHz"), "kHz")
    print(convert_time_to_freq(100, "ns", "MHz"), "MHz")

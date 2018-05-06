#!/usr/bin/env python


import numpy as np
from scipy.fftpack import fft


# rectify a signal for each point that goes beyond the specified (upper and
# lower) threshold
def rectify_signal(signal, up_thres=float("inf"), low_thres=0.0):
    for ind, val in enumerate(signal):
        if val < low_thres:
            signal[ind] = low_thres
        if val > up_thres:
            signal[ind] = up_thres
    return signal


# Fourier reconstruction / spectral analysis using l-BAM or l-AAM technique
# described in "A Poisson-spectral model for modelling temporal patterns in
# human data observed by a robot".
# @addition_method represents l-AAM technique
def reconstruct_signal(signal, addition_method=True):
    num_of_freqs = min(len(signal)/10, 15)
    if addition_method:
        spectrums, residue = get_accumulated_highest_n_freq(signal, num_of_freqs*2)
        spectrums = spectrums[0:num_of_freqs]
    else:
        residue = signal
        xf = np.linspace(0.0, len(signal), len(signal))
        spectrums = get_highest_n_freq(fft(signal), num_of_freqs)
        for [amp, phs, freq] in spectrums:
            wave = amp * np.cos((freq * 2.0 * np.pi * xf) + phs)
            residue -= wave
    reconstruction = 0
    for spectrum in spectrums:
        xf = np.linspace(0.0, len(signal), len(signal))
        wave = spectrum[0] * np.cos((spectrum[2]*2.0*np.pi*xf) + spectrum[1])
        reconstruction += wave
    return reconstruction, residue


# Addition Amplitude Model (l-AAM) technique to get the l highest frequencies.
# @max_addition is the maximum number of each frequency can be added up.
# @max_iteration is the maximum number of iteration to obtain the desired @num_of_freqs.
def get_accumulated_highest_n_freq(signal, num_of_freqs=15, max_addition=10, max_iteration=1000):
    xf = np.linspace(0.0, len(signal), len(signal))
    # initialise significant frequencies by taking frequency 0
    spectrums = fft(signal)
    [amp, phs, freq] = get_highest_n_freq(spectrums, 1)[0]
    frequencies = [[amp, phs, freq]]
    freq_counter = {freq: 1}
    exit_counter = 0
    while len(frequencies) < num_of_freqs:
        spectrums = fft(signal)
        # create a signal of the highest frequency
        # typically the highest frequency is frequency 0, that is why the second
        # highest is taken
        freqs = get_highest_n_freq(spectrums, 2)
        [amp, phs, freq] = freqs[1]
        if freq == 0:
            [amp, phs, freq] = freqs[0]
        wave = amp * np.cos((freq * 2.0 * np.pi * xf) + phs)
        # substracting signal with the wave
        signal -= wave
        if freq not in zip(*frequencies)[2]:
            frequencies.append([amp, phs, freq])
            freq_counter.update({freq: 1})
        else:
            for ind, val in enumerate(frequencies):
                if frequencies[ind][2] == freq and freq_counter[freq] < max_addition:
                    frequencies[ind][0] += amp
                    frequencies[ind][1] = ((
                        freq_counter[freq] * frequencies[ind][1]
                    ) + phs) / (freq_counter[freq] + 1)
                    freq_counter[freq] += 1
        exit_counter += 1
        if exit_counter >= max_iteration:
            break
    return frequencies, signal


# Best Amplitude Model (l-BAM) technique to get the l highest frequencies.
def get_highest_n_freq(freqs, n=15):
    N = len(freqs)
    freqs = freqs[0:N/2]
    indices = [i for i in range(len(freqs))]
    angles = np.angle(freqs)
    amplitudes = np.abs(freqs) / float(N)
    sorted_result = sorted(zip(amplitudes, angles, indices), reverse=True)
    n_freqs = sorted_result[:n]
    return n_freqs

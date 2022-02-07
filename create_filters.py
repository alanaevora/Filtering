import math
import numpy as np
from scipy import signal
import soundfile as sf
import sounddevice as sd


def white_noise(fs):
    # create white noise
    # taper white noise
    # pad white noise

    print("WHITE NOSIE\n")
    # create white noise
    wn = np.random.rand(2 * fs) - 0.5  # generate 2 sec of random noise (vector of evenly spaced random nums)
    
    # taper white noise
    window = np.ones(len(wn))  # array of 1s that is the length of the random noise
    ramp = np.linspace(0, np.pi, round(fs/5))   # creates a vector = linspace(min(x), max(x), #points)
    window[0:len(ramp)] = 0.5*np.sin(ramp-np.pi/2) + 0.5   # gradually taper at the beginning
    window[len(window) - len(ramp) : len(window)] = 0.5 * np.cos(ramp) + 0.5  # gradually taper at the end
    wn = wn * window

    # pad white noise
    padding = np.zeros(round(fs/6))  # array of 0s that is 1/6 of the length of the window
    pb_wn = []
    pb_wn = np.append(pb_wn, padding)
    pb_wn = np.append(pb_wn, wn)
    pb_wn = np.append(pb_wn, padding)
    pb_wn = np.append(pb_wn, padding)
    pb_wn = np.append(pb_wn, padding)
    # pb_wn = [padding, pb_wn, padding, padding, padding]
    sf.write('Noise_8ch.wav', pb_wn, fs)  # create wave file & save the random noise in it
    return pb_wn


def play_rec(fs, sound):
    # record to wave file source: https://python-sounddevice.readthedocs.io/en/0.4.1/usage.html#recording
    # possible source for APO: https://python-sounddevice.readthedocs.io/en/0.3.11/#sounddevice.AsioSetting
    file_reader, sample_rate = sf.read('Noise_8ch.wav')
    print("RECORDING...\n")
    r_wn = sd.playrec(file_reader, fs, channels=1)  # playing and recording simultaneously & saving audio as NumPy array
    sd.wait()  # wait for the recording to finish before moving on
    return r_wn

# def plot(spect, amp, name, x, y):
    # # Generate the figure **without using pyplot**.
    # print("GENERATING FIGURE FOR: " + name)
    # fig = Figure(figsize=(8, 4))
    # ax = fig.subplots()
    # ax.set_xlabel(x)
    # ax.set_ylabel(y)
    # ax.plot(spect, amp)
    # file = "static/images/"+ name + ".png"
    # fig.savefig(file)
    # return file


def calc_filter(iter, fs, nfft, pb_wn, r_wn):
    # adjust the amp
    # re-record with amp adjusted noise
    # equalize the length of recording and playback
    # calculate power spectrum
    # calculate the ratio of amplitudes at each frequency
    # use the ratio to design a filter
    # save the filter to a .wav file\=
    
    if iter == 1:
        pb_wn = pb_wn
    ##else:
        ## pb_wn = pb_filt_wn
    
    # adjust the amp
    if max(abs(r_wn)) > 0.2: # still need to figure out what this is doing
        print("-- exceeds max amp")
        adjust = 0.2 / max(abs(r_wn))
        pb_wn = adjust * pb_wn
    else:
        print("-- below max amp")
    
    # re-record with amp adjusted noise
    r_wn = play_rec(fs, pb_wn)

    # equalize the length of recording and playback
    if len(r_wn) > len(pb_wn):
        print("r_wn > pb_wn, shortening r_wn")
        r_wn[len(pb_wn)+1 : len(r_wn)] = []
    elif len(r_wn) < len(pb_wn):
        print("r_wn < pb_wn, lengthening r_wn (with silence)")
        r_wn = r_wn[len(r_wn):len(pb_wn)] + [0]*(len(pb_wn)-len(r_wn))

    print(len(r_wn), len(pb_wn))

    # calculate power spectrum
    pb_wn_spect, pb_wn_amps = signal.welch(pb_wn, nfft)
    r_wn_spect, r_wn_amps = signal.welch(r_wn, nfft) 

    # calculate the ratio of amplitudes at each frequency
    r = np.divide(pb_wn_amps, r_wn_amps) # the ratio between the original and the recorded

    print(pb_wn_spect)
    print(r_wn_spect)
    print(r)
    # use the ratio to design a filter
    # save the filter to a .wav file



if __name__ == "__main__":
    fs = 48000
    nfft = 1
    pb_wn = white_noise(fs)
    r_wn = play_rec(fs, pb_wn)

    iterations = [1]
    for iter in iterations:
        print("------ STARTING ITERATION #" + str(iter) + " ------")
        calc_filter(iter, fs, nfft, pb_wn, r_wn)
    
import numpy as np
import soundfile as sf
import sounddevice as sd
import matplotlib.pyplot as plt
from scipy import signal


def white_noise(fs):
    # create white noise
    # taper white noise
    # pad white noise

    print("WHITE NOSIE\n")
    # create white noise
    wn = np.random.rand(2 * fs) - 0.5  # generate 2 sec of random noise (vector of evenly spaced random nums)
    # plot_ta(wn, "wn: white noise")

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
    sf.write('white_noise.wav', pb_wn, fs)  #  save the random noise in .wav file
    return


def play_rec(rec_filename, play_filename):
    # record to wave file source: https://python-sounddevice.readthedocs.io/en/0.4.1/usage.html#recording
    # possible source for APO: https://python-sounddevice.readthedocs.io/en/0.3.11/#sounddevice.AsioSetting
    file_reader, fs = sf.read(play_filename)
    print("RECORDING...\n")
    recording = sd.playrec(file_reader, fs, channels=1)  # playing and recording simultaneously & saving audio as NumPy array
    sd.wait()  # wait for the recording to finish before moving on
    sf.write(rec_filename, recording, fs) #  save the recording in .wav file
    return


def plot_ta(x, title):
    # plots waveforms (time x amplitude)
    plt.plot(x)
    plt.title(title)
    plt.show()
    # the figure window blocks the rest of script from running until the window is clsoed
    # to fix this, you can run the plot function in a separate thread
    return

def plot_fa(f, psd, title):
    # plots ampltiude spectrums (frequency x power spectral density)
    plt.semilogy(f, psd) # semilogy transforms the y axis to log scaling: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.semilogy.html
    plt.xlabel('frequency [Hz]')
    plt.ylabel('amp')
    plt.title(title)
    plt.show()

def plot_fa3(f, psd1, psd2, psd3, title):
    # plots three ampltiude spectrums (frequency x power spectral density)
    plt.semilogy(f, psd1) 
    plt.semilogy(f, psd2)
    plt.semilogy(f, psd3)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('amp')
    plt.title(title)
    plt.show()

def plot_fa4(f, psd1, psd2, psd3, psd4, title):
    # plots three ampltiude spectrums (frequency x power spectral density)
    plt.semilogy(f, psd1) 
    plt.semilogy(f, psd2)
    plt.semilogy(f, psd3)
    plt.semilogy(f, psd4)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('amp')
    plt.title(title)
    plt.show()

def calc_filter(iter, fftSize):
    # adjust the amp of the playback
    # re-record with amp adjusted noise
    # equalize the length of recording and playback
    # calculate amplitude spectrum
    # calculate the ratio between the recording and the playback
    # use the ratio to design a filter
    # save the filter to a .wav file\=
    
    if iter == 1:
        recording_filename = "rec_white_noise" + str(iter) + ".wav"
        playback_filename = "white_noise.wav"
    # else:
        # recording_filename = "rec_white_noise" + iter + ".wav"
        # recording_filename = "rec_white_noise" + str(iter -1) + ".wav"
    
    play_rec(recording_filename, playback_filename)
    recording, r_fs = sf.read(recording_filename)
    playback, p_fs = sf.read(playback_filename)

    # adjust the amp
    if max(abs(recording)) > 0.2:
        print("-- exceeds max amp")
        adjust = 0.2 / max(abs(recording))
        playback = adjust * playback
        sf.write(playback_filename, playback, p_fs)
    else:
        print("-- below max amp")
    
    # re-record with amp adjusted noise
    play_rec(recording_filename, playback_filename)
    recording, r_fs = sf.read(recording_filename)
    playback, p_fs = sf.read(playback_filename)

    # equalize the length of recording and playback
    if len(recording) > len(playback):
        print("r_wn > pb_wn, shortening r_wn")
        recording[len(playback)+1 : len(recording)] = []
    elif len(recording) < len(playback):
        print("r_wn < pb_wn, lengthening r_wn (with silence)")
        recording = recording[len(recording):len(playback)] + [0]*(len(playback)-len(recording))

    # calculate ampltiude spectrums
    # signal.welch documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
    # freq_recording, Pxx_den_recording = signal.welch(recording, r_fs, nperseg=bin_length)
    # plot_fa(freq_recording, Pxx_den_recording, "psd recording")
    
    # freq_playback, Pxx_den_playback = signal.welch(playback, p_fs, nperseg=bin_length)
    # plot_fa(freq_playback, Pxx_den_playback, "psd playback")
    
    
    freq_recording, Pxx_den_recording = signal.welch(recording, r_fs, nfft=fftSize)
    plot_fa(freq_recording, Pxx_den_recording, "psd recording")
    
    freq_playback, Pxx_den_playback = signal.welch(playback, p_fs, nfft=fftSize)
    plot_fa(freq_playback, Pxx_den_playback, "psd playback")

    # calculate the ratio between the recording and the playback
    ratio = np.divide(Pxx_den_playback, Pxx_den_recording)
    plot_fa(freq_recording, ratio, "ratio")
    plot_fa3(freq_recording, Pxx_den_recording, Pxx_den_playback, ratio, "rec, playback, ratio")

    # use the ratio to design a filter
    filter_design = signal.firwin2(fftSize + 1, freq_playback, ratio, fs=fs)
    freq_design, Pxx_den_design = signal.welch(filter_design, r_fs, nfft=fftSize)
    plot_fa3(freq_design, Pxx_den_recording, Pxx_den_design, ratio, "rec, filter design, ratio")

    container = np.zeros(len(filter_design))
    container[0] = 1 # look at the dimensions

    filtered_data = signal.lfilter(filter_design, container, recording)
    
    freq_filtered, Pxx_den_filtered = signal.welch(filtered_data, r_fs, nfft=fftSize)
    plot_fa(freq_filtered, Pxx_den_filtered, "psd filtered")

    plot_fa4(freq_recording, Pxx_den_recording, Pxx_den_playback,ratio, Pxx_den_filtered, "all four")
    # save the filter to a .wav file
    

if __name__ == "__main__":
    fs = 48000
    fftSize = 16768 # fftLength in Rex's code, nfft in signal.welch documentation
    white_noise(fs)

    iterations = [1]
    for iter in iterations:
        print("------ STARTING ITERATION #" + str(iter) + " ------")
        calc_filter(iter, fftSize)

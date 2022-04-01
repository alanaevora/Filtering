"""
This is a Python script that generates compensation filters.
"""

import numpy as np
import soundfile as sf
import sounddevice as sd
import matplotlib.pyplot as plt
from scipy import signal

"""******************* BELOW ARE "HELPER" PLOT FUNCTIONS. THEY CAN BE USED TO PLOT AUDIO FOR TROUBLESHOOTING. *******************"""

"""
This is a function that plots waveforms together (if there are multiple), using time and amplitude domains. 
    Name: plot_ta
    Input:  sounds_list (sounds to plot, as a list)
            title (title of plot, as a string)
    Output: none

    Call Example: time_amp_plots([sound1, sound2, sound3], "This is a title.")
"""
def time_amp_plots(sounds_list, title):
    for sound in sounds_list:
        plt.plot(sound)
    plt.title(title)
    plt.show()

"""
This is a function that plots amplitude spectrums on the same window, using frequency and amplitude domains. 
    Name: freq_amp_plots
    Input:  freq_points (frequency points, as a numpy array) (all of the sounds must have this in common)
            amp_points_list (list of amp_points, as a list)
            title (title of plot, as a string)
    Output: none

    Call Example: freq_amp_plots(frequency_points, [amp_points1, amp_points2, amp_points3], "This is a title.")
"""
def freq_amp_plots(freq_points, amp_points_list, title):
    for amp_points in amp_points_list:
        plt.semilogy(freq_points, amp_points)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('amplitude [dB]')
    plt.title(title)
    plt.show()


"""******************* BELOW ARE COMPUTING FUNCTIONS. THIS IS THE BULK OF THE SCRIPT. *******************"""

"""
This function creates a wav file containing white noise (white_noise.wav).
    Name:   white_noise
    Input:  fs (the sampling rate, as an int)
    Output: none

    Call Example: white_noise("41000")

    STEPS:  1. create white noise
            2. taper white noise
            3. pad white noise

    SOURCES:
    MATLAB SCRIPT SOURCE
"""
def white_noise(fs):
    print("\n***  white_noise ***")
    print("1. create white noise")
    wn = np.random.rand(2 * fs) - 0.5  # creates a vector of evenly spaced random numbers (aka 2 seconds of white noise)

    print("2. taper white noise")
    window = np.ones(len(wn))  # creates a vector of 1s that is the same length of the random noise
    ramp = np.linspace(0, np.pi, round(fs/5))   # creates a vector later used modify the window to taper the white noise | np.linspace(min(x), max(x), #points)
    window[0:len(ramp)] = 0.5*np.sin(ramp-np.pi/2) + 0.5   # gradually taper window at the beginning
    window[len(window) - len(ramp) : len(window)] = 0.5 * np.cos(ramp) + 0.5  # gradually taper window at the end
    wn = wn * window # apply the taper window to the white noise

    print("pad white noise")
    padding = np.zeros(round(fs/6))  # vector of 0s that is 1/6 of the length of the window
    pb_wn = []
    pb_wn = np.append(pb_wn, padding)
    pb_wn = np.append(pb_wn, wn)
    pb_wn = np.append(pb_wn, padding)
    pb_wn = np.append(pb_wn, padding)
    pb_wn = np.append(pb_wn, padding)
    time_amp_plots([pb_wn], "pbwn") # <- uncomment this if you want to plot the white noise
    sf.write('white_noise.wav', pb_wn, fs)  # save the random noise in .wav file

"""
This is a function that plays & records audio simultaneously and saves the recording to a wav file. 
    Name:   play_rec
    Input:  rec_filename (name of file to create for the recording, as a string)
            play_filename (name of file to play, as a string)
    Output: none

    Call Example: play_rec("recording.wav", "play.wav")

    SOURCES:
    simultaneous playback and recording: https://python-sounddevice.readthedocs.io/en/0.4.1/usage.html#simultaneous-playback-and-recording
"""
def play_rec(rec_filename, play_filename):
    print("<<recording>>")
    file_reader, fs = sf.read(play_filename)
    recording = sd.playrec(file_reader, fs, channels=1)  # playing and recording simultaneously & saving audio as numpy array
    sd.wait()  # wait for the recording to finish before moving on
    sf.write(rec_filename, recording, fs) #  save the recording in .wav file

"""
This is a function that calculates compensation filters and saves them as .wav file.  
    Name:   calc_filter
    Input:  iter (the iteration number, as an int)
            fftSize (fast fourier transform size, as an int)
    Output: none

    Call Example: calc_filter(3, 16768)

    STEPS:
    1. adjust the amp of the playback
    2. re-record with amp adjusted noise
    3. equalize the length of recording and playback
    4. calculate amplitude spectrum
    5. calculate the ratio between the recording and the playback
    6. use the ratio to design a filter
    7. save the filter to a .wav file

    SOURCES:
    signal.welch (generating amplitude spectrums): https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
"""


def calc_filter(iter, fftSize):
    print("***  calc_filters    ***")
    if iter == 1:
        recording_filename = "rec_white_noise" + str(iter) + ".wav"
        playback_filename = "white_noise.wav"
    # else:
        # recording_filename = "rec_white_noise" + iter + ".wav"
        # recording_filename = "rec_white_noise" + str(iter -1) + ".wav"
    
    play_rec(recording_filename, playback_filename)
    recording, r_fs = sf.read(recording_filename)
    playback, p_fs = sf.read(playback_filename)

    print("1. adjust the amp of the playback")
    if max(abs(recording)) > 0.2:
        print("-- exceeds max amp")
        adjust = 0.2 / max(abs(recording))
        playback = adjust * playback
        sf.write(playback_filename, playback, p_fs)
    else:
        print("-- below max amp")
    
    print("2. re-record with amp adjusted noise")
    play_rec(recording_filename, playback_filename)
    recording, r_fs = sf.read(recording_filename)
    playback, p_fs = sf.read(playback_filename)

    print("3. equalize the length of recording and playback")
    if len(recording) > len(playback):
        print("r_wn > pb_wn, shortening r_wn")
        recording[len(playback)+1 : len(recording)] = []
    elif len(recording) < len(playback):
        print("r_wn < pb_wn, lengthening r_wn (with silence)")
        recording = recording[len(recording):len(playback)] + [0]*(len(playback)-len(recording))

    print("4. calculate ampltiude spectrums")    
    freq_points_rec, amp_points_rec = signal.welch(recording, r_fs, nfft=fftSize)
    freq_points_playback, amp_points_playback = signal.welch(playback, p_fs, nfft=fftSize)
    # freq_amp_plots(freq_points_rec, [amp_points_rec], "amp spec recording") # <- uncomment this if you want to plot the amplitude spectrum of the recording
    # freq_amp_plots(freq_points_playback, [amp_points_playback], "amp spec playback")  # <- uncomment this if you want to plot the amplitude spectrum of the playback

    print("5. calculate the ratio between the recording and the playback")
    amp_points_ratio = np.divide(amp_points_playback, amp_points_rec)
    # freq_amp_plots(freq_points_rec, [amp_points_ratio], "ratio") #<- uncomment this if you want to plot the amplitude spectrum of the ratio
    # freq_amp_plots(freq_points_rec, [amp_points_rec, amp_points_playback, amp_points_ratio], "rec, playback, ratio") # <- uncomment this if you want to plot the recording, playback, and ratio amplitude spectrums together

    print("6. use the ratio to design a filter")
    # filter_design = signal.firwin2(fftSize + 1, freq_playback, ratio, fs=fs)
    # freq_design, Pxx_den_design = signal.welch(filter_design, r_fs, nfft=fftSize)
    # plot_fa3(freq_design, Pxx_den_recording, Pxx_den_design, ratio, "rec (b), filter design (o), ratio(g)")

    # container = np.zeros(len(filter_design))
    # container[0] = 1 # look at the dimensions

    # filtered_data = signal.lfilter(filter_design, container, recording)
    
    # freq_filtered, Pxx_den_filtered = signal.welch(filtered_data, r_fs, nfft=fftSize)
    # plot_fa(freq_filtered, Pxx_den_filtered, "psd filtered")

    # plot_fa4(freq_recording, Pxx_den_recording, Pxx_den_playback, ratio, Pxx_den_filtered, "rec (b), playback (o), ratio (g), filtered (r)")
    # # save the filter to a .wav file


"""******************* BELOW IS THE MAIN THREAD. THIS IS WHAT ACTUALLY RUNS WHEN YOU RUN THIS PYTHON FILE. *******************"""

if __name__ == "__main__":
    fs = 48000 # the sampling frequency
    fft_size = 16768 # the fast fourier transform size. this is fftLength in Rex's code
    white_noise(fs) #  calls the white_noise function

    iterations = [1] # the number of iterations to perfrom (for testing right now, it is 1 iteraction. usually it is 4 iterations)
    for iter in iterations:
        print("\n------ STARTING ITERATION #" + str(iter) + " ------")
        calc_filter(iter, fft_size) # calls the calc_filter function

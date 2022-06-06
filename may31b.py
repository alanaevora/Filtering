import numpy as np
import soundfile as sf
import sounddevice as sd
import matplotlib.pyplot as plt
from scipy import signal
import sys

fs = 48000
iter = 1
fft = 256
lo_val = 100
hi_val = 20000

print("1. create white noise")

print("2. determine iteration / filter cycle")

print("3. play and record")
recording, r_fs = sf.read("rec_white_noise1.wav")   # gets the recording
playback, p_fs = sf.read("white_noise.wav")         # get the playback

print("4. equalize the length of recording and playback")
print("pre len rec: ", len(recording))
print("pre len pb: ", len(playback))
if len(recording) > len(playback):                          # if the recording is longer than the playback, shorten the recording
    print("rec > pb")
    recording[len(playback)+1 : len(recording)] = []
elif len(recording) < len(playback):                        # if the playback is longer than the recording, lengthen the recording with silence
    print("rec < pb")
    recording = recording[len(recording):len(playback)] + [0]*(len(playback)-len(recording))

print("post len rec: ", len(recording))
print("post len pb: ", len(playback))

print("5. plot waveform of recorded noise")

print("6.get amp spectrums")
freq_rec, amp_rec = signal.welch(recording, r_fs, window='hamming', nperseg=fft, scaling='spectrum', detrend=False)
freq_pb, amp_pb = signal.welch(playback, p_fs, window='hamming', nperseg=fft, scaling='spectrum', detrend=False)
amp_rec = np.power(amp_rec, .5)
amp_pb = np.power(amp_pb, .5)

np.set_printoptions(threshold=sys.maxsize)
print(amp_pb)
print(amp_rec)

print("7. calculate the ratio between rec and pb")
lo = int(np.floor(lo_val / (fs/fft)))
hi = int(np.floor(hi_val / (fs/fft)))
amp_ratio = np.divide(amp_pb, amp_rec)
amp_ratio[0:lo] = 0
amp_ratio[hi+10:len(amp_ratio)] = 0

print("RRRRRRR")
print(amp_ratio)
plt.semilogy(freq_pb, amp_pb)
plt.semilogy(freq_pb, amp_rec)
plt.plot(freq_pb, amp_ratio)
plt.show()

print("8. design the filter")
amp_ratio[0] = 0 # maybe matlab does this behind the sences to correct?
freq = np.linspace(0,1,len(amp_ratio))
freq = np.multiply(fs/2, freq)

b_coefficient = signal.firwin2(fft, freq, gain_ratio, fs=fs)

print("BBBBBBBBBB")
print(b_coefficient)

freq_b, amp_b= signal.welch(b_coefficient, p_fs, window='hamming', nperseg=fft, scaling='spectrum', detrend=False)
amp_filter = np.power(amp_b, .5)
plt.semilogy(freq_pb, amp_ratio)
plt.semilogy(freq_b, amp_b)
plt.show()

print("9. filter the original (or previously filtered) noise")
a_coefficient = np.zeros(len(b_coefficient))
a_coefficient[0] = 1 # look at the dimensions

# plt.semilogy(signal.freqz(b_coefficient, a_coefficient))
# plt.show()

filtered_signal = signal.lfilter(b_coefficient, a_coefficient, recording)

freq_filtered, amp_filtered = signal.welch(filtered_signal, r_fs, nfft=fft)

freq_pb = freq_pb[1:hi]
amp_ratio = amp_ratio[1:hi]
amp_pb = amp_pb[1:hi]
amp_filtered = amp_filtered[1:hi]
plt.semilogy(freq_pb, amp_ratio)
plt.semilogy(freq_pb, amp_filtered)
plt.semilogy(freq_pb, amp_pb)
plt.show()
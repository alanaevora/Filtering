import soundfile as sf
import numpy as np
import matplotlib.pyplot as plt 

pb_peak_amp = 0.5
target_amp = 0.8
seglen = 8192   # length of segments that are used for amp adjustments
fs = 48000

print("1. get filters and stimulus files")
print("2. filter stimulus")
pb_filt_stim, fs_pb_filt_stim = sf.read("pb_filt_stim.wav")

print("\n3. play & record filtered stimulus")
rec_filt_stim, fs_rec_filt_stim = sf.read("rec_filt_stim.wav")   # get the filtered stimulus recording

print("\n4. get time delay") # PYTHON V MATLAB TESTED
peak_rec = max(abs(rec_filt_stim))                      # the max amp value in recording
peak_loc = int(np.where(rec_filt_stim == peak_rec)[0][0])    # finds the first x (time) location of the maximum y (amp) value
pb_vs_rec = np.correlate(pb_filt_stim, rec_filt_stim, "full")
peak_corr = int(np.where(abs(pb_vs_rec) == max(abs(pb_vs_rec)))[0][0])
time_delay = (len(pb_vs_rec)/2) - peak_corr
peak_loc = peak_loc - time_delay

print("peak_corr: ", peak_corr)
print("time_delay: ", time_delay)
print("peak_loc: ", peak_loc)

print("\n5. get and window segments of recordings for calibration")
seg = pb_filt_stim[round(peak_loc - seglen/2) + 1 : round(peak_loc + seglen/2)] # extract segment (points before peak, peak, points after peak)
window = np.ones(len(seg)) # create a vector of 1s that is the same length as the segments
ramp = np.linspace(0, np.pi, round(seglen/8)) # create a vector later used to modify the window to taper the segments
window[0:len(ramp)] = 0.5*np.sin(ramp-np.pi/2) + 0.5 # gradually taper window at the beginning 
window[len(window) - len(ramp) : len(window)] = 0.5*np.cos(ramp) + 0.5 # gradually taper window at the end
seg = seg * window # apply the window
seg_max = max(abs(seg)) # get the max amp of windowed segments

# plt.plot(seg)
plt.show()
print("sigMax: ", seg_max)

print("\n6. generate ramped segments of increasing amp")
step_max = 1.25 * seg_max / (peak_rec / pb_peak_amp)    # ensures that the steps go above the max amp
steps = np.linspace(0.01, step_max, 20)                 # create a vector 
seg_ramp = []
pb_seg_locs = []
step_counter = 0
for i in steps:
    seg_ramp = np.append(seg_ramp, np.multiply(steps[step_counter],seg))  # multiplies each step by the segment and adds it to the array, creating an increasing ramp of the segment
    pb_seg_locs.append([step_counter*seglen + 1, (step_counter + 1)*seglen]) # saves the start x value and end x value of each ramp in the ramp [(start, end), (start, end) ...] with the length of the padding already accounted for
    step_counter = step_counter + 1

print(steps)

pb_seg_ramp = []
padding = np.zeros(round(fs/2))
pb_seg_ramp = np.append(pb_seg_ramp, padding)
pb_seg_ramp = np.append(pb_seg_ramp, seg_ramp)
pb_seg_ramp = np.append(pb_seg_ramp, padding)
pb_seg_ramp = np.append(pb_seg_ramp, padding)
pb_seg_ramp = np.append(pb_seg_ramp, padding)     # add padding (silence) before and after the ramp

for x in pb_seg_locs:
    x[0] = x[0] + len(padding)
    x[1] = x[1] + len(padding)

print("sig_locs: ")
print(pb_seg_locs)

print("\n7. play & record ramped segments")
rec_seg_ramp, fs_rec_seg_ramp = sf.read("rec_seg_ramp.wav") # save the recording ramp in .wav file

print("\n8. get time delay") # CHECK PYTHON VS MATLAB
rec_ramp_peak = max(abs(rec_seg_ramp)) # max amp value
rec_ramp_peak_loc = int(np.where(abs(rec_seg_ramp) == rec_ramp_peak)[0][0]) # location of max amp value in ramp rec (middle of the last segment)
shift = round(rec_ramp_peak_loc + seglen/2) - pb_seg_locs[19][1] # adds 1/2 length of segment to get to end of the ramp, subtracts the location of the end of the ramp in pb to find difference between rec & pb (which is equivalent to the time delay)

rec_seg_locs = []
for y in pb_seg_locs:
    start = y[0] + shift
    end = y[1] + shift
    rec_seg_locs.append([start, end])

print("shift: ", shift)
print("rec_sig_locs: ")
print(rec_seg_locs)

print("\n9. calculate playback amp needed to reach target amp")
pks = []
for j in rec_seg_locs:
    thisseg = rec_seg_ramp[j[0]:j[1]] # gets each segment from recording
    pks.append(max(abs(thisseg))) # gets max amp value, and adds to list

p = np.polyfit(pks, steps, 1) # returns slope (degree of 0) and y intercept (degree of 1) of line of best fit
mult = np.polyval(p, target_amp) # plugs target amp in for x to find y (the amp to play at to get the target amp)

plt.plot(pks,steps)
plt.show()

print(p)

print("pks: ")
print(pks)
print("p: ", p)
print("mult ", mult)
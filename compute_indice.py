#!/usr/bin/env python

"""
    Set of functions to compute acoustic indices in the framework of Soundscape Ecology.

    Some features are inspired or ported from those proposed in:
        - seewave R package (http://rug.mnhn.fr/seewave/) / Jerome Sueur, Thierry Aubin and  Caroline Simonis
        - soundecology R package (http://cran.r-project.org/web/packages/soundecology/index.html) / Luis J. Villanueva-Rivera and Bryan C. Pijanowski

    This file use an object oriented type for audio files described in the file "acoustic_index.py".



"""

__author__ = "Patrice Guyot"
__version__ = "0.4"
__credits__ = ["Patrice Guyot", "Alice Eldridge", "Mika Peck"]
__email__ = ["guyot.patrice@gmail.com", "alicee@sussex.ac.uk", "m.r.peck@sussex.ac.uk"]
__status__ = "Development"


from scipy import signal, fftpack
import numpy as np
import matplotlib.pyplot as plt



#--------------------------------------------------------------------------------------------------------------------------------------------------
def compute_spectrogram(file, windowLength=512, windowHop= 256, scale_audio=True, square=True, windowType='hanning', centered=False, normalized = False ):
    """
    Compute a spectrogram of an audio signal.
    Return a list of list of values as the spectrogram, and a list of frequencies.

    Keyword arguments:
    file -- the real part (default 0.0)

    Parameters:
    file: an instance of the AudioFile class.
    windowLength: length of the fft window (in samples)
    windowHop: hop size of the fft window (in samples)
    scale_audio: if set as True, the signal samples are scale between -1 and 1 (as the audio convention). If false the signal samples remains Integers (as output from scipy.io.wavfile)
    square: if set as True, the spectrogram is computed as the square of the magnitude of the fft. If not, it is the magnitude of the fft.
    hamming: if set as True, the spectrogram use a correlation with a hamming window.
    centered: if set as true, each resulting fft is centered on the corresponding sliding window
    normalized: if set as true, divide all values by the maximum value
    """

    if scale_audio:
        sig = file.sig_float # use signal with float between -1 and 1
    else:
        sig = file.sig_int # use signal with integers

    W = signal.get_window(windowType, windowLength, fftbins=False)
    halfWindowLength = int(windowLength/2)

    if centered:
        time_shift = int(windowLength/2)
        times = range(time_shift, len(sig)+1-time_shift, windowHop) # centered
        frames = [sig[i-time_shift:i+time_shift]*W for i in times] # centered frames
    else:
        times = range(0, len(sig)-windowLength+1, windowHop)
        frames = [sig[i:i+windowLength]*W for i in times]

    if square:
        spectro =  [abs(np.fft.rfft(frame, windowLength))[0:halfWindowLength]**2 for frame in frames]
    else:
        spectro =  [abs(np.fft.rfft(frame, windowLength))[0:halfWindowLength] for frame in frames]

    spectro=np.transpose(spectro) # set the spectro in a friendly way

    if normalized:
        spectro = spectro/np.max(spectro) # set the maximum value to 1 y

    frequencies = [e * file.niquist / float(windowLength / 2) for e in range(halfWindowLength)] # vector of frequency<-bin in the spectrogram
    return spectro, frequencies





#-----------------------------------------------------------------------------
def compute_ACI(spectro,j_bin):
    """
    Compute the Acoustic Complexity Index from the spectrogram of an audio signal.

    Reference: Pieretti N, Farina A, Morri FD (2011) A new methodology to infer the singing activity of an avian community: the Acoustic Complexity Index (ACI). Ecological Indicators, 11, 868-873.

    Ported from the soundecology R package.

    spectro: the spectrogram of the audio signal
    j_bin: temporal size of the frame (in samples)


    """

    #times = range(0, spectro.shape[1], j_bin) # relevant time indices
    times = range(0, spectro.shape[1]-10, j_bin) # alternative time indices to follow the R code

    jspecs = [np.array(spectro[:,i:i+j_bin]) for i in times]  # sub-spectros of temporal size j

    aci = [sum((np.sum(abs(np.diff(jspec)), axis=1) / np.sum(jspec, axis=1))) for jspec in jspecs] 	# list of ACI values on each jspecs
    main_value = sum(aci)
    temporal_values = aci

    return main_value, temporal_values # return main (global) value, temporal values





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_BI(spectro, frequencies, min_freq = 2000, max_freq = 8000):
    """
    Compute the Bioacoustic Index from the spectrogram of an audio signal.
    In this code, the Bioacoustic Index correspond to the area under the mean spectre (in dB) minus the minimum frequency value of this mean spectre.

    Reference: Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne remote sensing. Ecological Applications 17: 2137-2144.

    spectro: the spectrogram of the audio signal
    frequencies: list of the frequencies of the spectrogram
    min_freq: minimum frequency (in Hertz)
    max_freq: maximum frequency (in Hertz)

    Ported from the soundecology R package.
    """

    min_freq_bin = int(np.argmin([abs(e - min_freq) for e in frequencies])) # min freq in samples (or bin)
    max_freq_bin = int(np.ceil(np.argmin([abs(e - max_freq) for e in frequencies]))) # max freq in samples (or bin)

    min_freq_bin = min_freq_bin - 1 # alternative value to follow the R code



    spectro_BI = 20 * np.log10(spectro/np.max(spectro))  #  Use of decibel values. Equivalent in the R code to: spec_left <- spectro(left, f = samplingrate, wl = fft_w, plot = FALSE, dB = "max0")$amp
    spectre_BI_mean = 10 * np.log10 (np.mean(10 ** (spectro_BI/10), axis=1))     # Compute the mean for each frequency (the output is a spectre). This is not exactly the mean, but it is equivalent to the R code to: return(a*log10(mean(10^(x/a))))
    spectre_BI_mean_segment =  spectre_BI_mean[min_freq_bin:max_freq_bin]   # Segment between min_freq and max_freq
    spectre_BI_mean_segment_normalized = spectre_BI_mean_segment - min(spectre_BI_mean_segment) # Normalization: set the minimum value of the frequencies to zero.
    area = np.sum(spectre_BI_mean_segment_normalized / (frequencies[1]-frequencies[0]))   # Compute the area under the spectre curve. Equivalent in the R code to: left_area <- sum(specA_left_segment_normalized * rows_width)

    return area

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_SH(spectro):
    """
    Compute Spectral Entropy of Shannon from the spectrogram of an audio signal.

    spectro: the spectrogram of the audio signal

    Ported from the seewave R package.
    """
    N = spectro.shape[0]
    spec = np.sum(spectro,axis=1)
    spec = spec / np.sum(spec)  # Normalization by the sum of the values
    main_value = - sum([y * np.log2(y) for y in spec]) / np.log2(N)  #Equivalent in the R code to: z <- -sum(spec*log(spec))/log(N)
    #temporal_values = [- sum([y * np.log2(y) for y in frame]) / (np.sum(frame) * np.log2(N)) for frame in spectro.T]
    return main_value





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_TH(file, integer=True):
    """
    Compute Temporal Entropy of Shannon from an audio signal.

    file: an instance of the AudioFile class.
    integer: if set as True, the Temporal Entropy will be compute on the Integer values of the signal. If not, the signal will be set between -1 and 1.

    Ported from the seewave R package.
    """
    if integer:
        sig=file.sig_int
    else:
        sig=file.sig_float

    #env = abs(signal.hilbert(sig)) # Modulo of the Hilbert Envelope
    env = abs(signal.hilbert(sig, fftpack.helper.next_fast_len(len(sig)))) # Modulo of the Hilbert Envelope, computed with the next fast length window

    env = env / np.sum(env)  # Normalization
    N = len(env)
    return - sum([y * np.log2(y) for y in env]) / np.log2(N)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_NDSI(file, windowLength = 1024, anthrophony=[1000,2000], biophony=[2000,11000]):
    """
    Compute Normalized Difference Sound Index from an audio signal.
    This function compute an estimate power spectral density using Welch's method.

    Reference: Kasten, Eric P., Stuart H. Gage, Jordan Fox, and Wooyeong Joo. 2012. The Remote Environ- mental Assessment Laboratory's Acoustic Library: An Archive for Studying Soundscape Ecology. Ecological Informatics 12: 50-67.

    windowLength: the length of the window for the Welch's method.
    anthrophony: list of two values containing the minimum and maximum frequencies (in Hertz) for antrophony.
    biophony: list of two values containing the minimum and maximum frequencies (in Hertz) for biophony.

    Inspired by the seewave R package, the soundecology R package and the original matlab code from the authors.
    """

    #frequencies, pxx = signal.welch(file.sig_float, fs=file.sr, window='hamming', nperseg=windowLength, noverlap=windowLength/2, nfft=windowLength, detrend=False, return_onesided=True, scaling='density', axis=-1) # Estimate power spectral density using Welch's method
    # TODO change of detrend for apollo
    frequencies, pxx = signal.welch(file.sig_float, fs=file.sr, window='hamming', nperseg=windowLength, noverlap=windowLength/2, nfft=windowLength, detrend='constant', return_onesided=True, scaling='density', axis=-1) # Estimate power spectral density using Welch's method
    avgpow = pxx * frequencies[1] # use a rectangle approximation of the integral of the signal's power spectral density (PSD)
    #avgpow = avgpow / np.linalg.norm(avgpow, ord=2) # Normalization (doesn't change the NDSI values. Slightly differ from the matlab code).

    min_anthro_bin=np.argmin([abs(e - anthrophony[0]) for e in frequencies])  # min freq of anthrophony in samples (or bin) (closest bin)
    max_anthro_bin=np.argmin([abs(e - anthrophony[1]) for e in frequencies])  # max freq of anthrophony in samples (or bin)
    min_bio_bin=np.argmin([abs(e - biophony[0]) for e in frequencies])  # min freq of biophony in samples (or bin)
    max_bio_bin=np.argmin([abs(e - biophony[1]) for e in frequencies])  # max freq of biophony in samples (or bin)

    anthro = np.sum(avgpow[min_anthro_bin:max_anthro_bin])
    bio = np.sum(avgpow[min_bio_bin:max_bio_bin])

    ndsi = (bio - anthro) / (bio + anthro)
    return ndsi


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def gini(values):
    """
    Compute the Gini index of values.

    values: a list of values

    Inspired by http://mathworld.wolfram.com/GiniCoefficient.html and http://en.wikipedia.org/wiki/Gini_coefficient
    """
    y = sorted(values)
    n = len(y)
    G = np.sum([i*j for i,j in zip(y,range(1,n+1))])
    G = 2 * G / np.sum(y) - (n+1)
    return G/n


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_AEI(spectro, freq_band_Hz, max_freq=10000, db_threshold=-50, freq_step=1000):
    """
    Compute Acoustic Evenness Index of an audio signal.

    Reference: Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246.

    spectro: spectrogram of the audio signal
    freq_band_Hz: frequency band size of one bin of the spectrogram (in Hertz)
    max_freq: the maximum frequency to consider to compute AEI (in Hertz)
    db_threshold: the minimum dB value to consider for the bins of the spectrogram
    freq_step: size of frequency bands to compute AEI (in Hertz)

    Ported from the soundecology R package.
    """

    bands_Hz = range(0, max_freq, freq_step)
    bands_bin = [f / freq_band_Hz for f in bands_Hz]

    spec_AEI = 20*np.log10(spectro/np.max(spectro))
    spec_AEI_bands = [spec_AEI[int(bands_bin[k]):int(bands_bin[k]+bands_bin[1]),] for k in range(len(bands_bin))]

    values = [np.sum(spec_AEI_bands[k]>db_threshold)/float(spec_AEI_bands[k].size) for k in range(len(bands_bin))]

    return gini(values)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_ADI(spectro, freq_band_Hz,  max_freq=10000, db_threshold=-50, freq_step=1000):
    """
    Compute Acoustic Diversity Index.

    Reference: Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246.

    spectro: spectrogram of the audio signal
    freq_band_Hz: frequency band size of one bin of the spectrogram (in Hertz)
    max_freq: the maximum frequency to consider to compute ADI (in Hertz)
    db_threshold: the minimum dB value to consider for the bins of the spectrogram
    freq_step: size of frequency bands to compute ADI (in Hertz)


    Ported from the soundecology R package.
    """


    bands_Hz = range(0, max_freq, freq_step)
    bands_bin = [f / freq_band_Hz for f in bands_Hz]

    spec_ADI = 20*np.log10(spectro/np.max(spectro))
    spec_ADI_bands = [spec_ADI[int(bands_bin[k]):int(bands_bin[k]+bands_bin[1]),] for k in range(len(bands_bin))]

    values = [np.sum(spec_ADI_bands[k]>db_threshold)/float(spec_ADI_bands[k].size) for k in range(len(bands_bin))]

    # Shannon Entropy of the values
    #shannon = - sum([y * np.log(y) for y in values]) / len(values)  # Follows the R code. But log is generally log2 for Shannon entropy. Equivalent to shannon = False in soundecology.

    # The following is equivalent to shannon = True (default) in soundecology. Compute the Shannon diversity index from the R function diversity {vegan}.
    #v = [x/np.sum(values) for x in values]
    #v2 = [-i * j  for i,j in zip(v, np.log(v))]
    #return np.sum(v2)

    # Remove zero values (Jan 2016)
    values = [value for value in values if value != 0]

    #replace zero values by 1e-07 (closer to R code, but results quite similars)
    #values = [x if x != 0 else 1e-07 for x in values]

    return np.sum([-i/ np.sum(values) * np.log(i / np.sum(values))  for i in values])




#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_zcr(file, windowLength=512, windowHop= 256):
    """
    Compute the Zero Crossing Rate of an audio signal.

    file: an instance of the AudioFile class.
    windowLength: size of the sliding window (samples)
    windowHop: size of the lag window (samples)

    return: a list of values (number of zero crossing for each window)
    """


    sig = file.sig_int # Signal on integer values

    times = range(0, len(sig)- windowLength +1, windowHop)
    frames = [sig[i:i+windowLength] for i in times]
    return [len(np.where(np.diff(np.signbit(x)))[0])/float(windowLength) for x in frames]


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_rms_energy(file, windowLength=512, windowHop=256, integer=False):
    """
    Compute the RMS short time energy.

    file: an instance of the AudioFile class.
    windowLength: size of the sliding window (samples)
    windowHop: size of the lag window (samples)
    integer: if set as True, the Temporal Entropy will be compute on the Integer values of the signal. If not, the signal will be set between -1 and 1.

    return: a list of values (rms energy for each window)
    """
    if integer:
        sig=file.sig_int
    else:
        sig=file.sig_float

    times = range(0, len(sig) - windowLength+1, windowHop)
    frames = [sig[i:i + windowLength] for i in times]
    return [np.sqrt(sum([ x**2 for x in frame ]) / windowLength) for frame in frames]


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_spectral_centroid(spectro, frequencies):
    """
    Compute the spectral centroid of an audio signal from its spectrogram.

    spectro: spectrogram of the audio signal
    frequencies: list of the frequencies of the spectrogram
    """

    centroid = [ np.sum(magnitudes*frequencies) / np.sum(magnitudes) for magnitudes in spectro.T]


    return centroid



#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_wave_SNR(file, frame_length_e=512, min_DB=-60, window_smoothing_e=5, activity_threshold_dB=3, hist_number_bins = 100, dB_range = 10, N = 0):
    """

    Computes indices from the Signal to Noise Ratio of a waveform.

    file: an instance of the AudioFile class.
    window_smoothing_e: odd number for sliding mean smoothing of the histogram (can be 3, 5 or 7)
    hist_number_bins - Number of columns in the histogram
    dB_range - dB range to consider in the histogram
    N: The decibel threshold for the waveform is given by the modal intensity plus N times the standard deviation. Higher values of N will remove more energy from the waveform.

    Output:
        Signal-to-noise ratio (SNR): the decibel difference between the maximum envelope amplitude in any minute segment and the background noise.
        Acoustic activity: the fraction of frames within a one minute segment where the signal envelope is more than 3 dB above the level of background noise
        Count of acoustic events: the number of times that the signal envelope crosses the 3 dB threshold
        Average duration of acoustic events: an acoustic event is a portion of recordingwhich startswhen the signal envelope crosses above the 3 dB threshold and ends when it crosses belowthe 3 dB threshold.

    Ref: Towsey, Michael W. (2013) Noise removal from wave-forms and spectro- grams derived from natural recordings of the environment.
    Towsey, Michael (2013), Noise Removal from Waveforms and Spectrograms Derived from Natural Recordings of the Environment. Queensland University of Technology, Brisbane.
    """

    half_window_smoothing = int(window_smoothing_e/2)

    times = range(0, len(file.sig_int)-frame_length_e+1, frame_length_e)
    wave_env = 20*np.log10([np.max(abs(file.sig_float[i : i + frame_length_e])) for i in times])

    minimum = np.max((np.min(wave_env), min_DB)) # If the minimum value is less than -60dB, the minimum is set to -60dB

    hist, bin_edges = np.histogram(wave_env, range=(minimum, minimum + dB_range), bins=hist_number_bins, density=False)


    hist_smooth = ([np.mean(hist[i - half_window_smoothing: i + half_window_smoothing]) for i in range(half_window_smoothing, len(hist) - half_window_smoothing)])
    hist_smooth = np.concatenate((np.zeros(half_window_smoothing), hist_smooth, np.zeros(half_window_smoothing)))

    modal_intensity = np.argmax(hist_smooth)

    if N>0:
        count_thresh = 68 * sum(hist_smooth) / 100
        count = hist_smooth[modal_intensity]
        index_bin = 1
        while count < count_thresh:
            if modal_intensity + index_bin <= len(hist_smooth):
                count = count + hist_smooth[modal_intensity + index_bin]
            if modal_intensity - index_bin >= 0:
                count = count + hist_smooth[modal_intensity - index_bin]
            index_bin += 1
        thresh = np.min((hist_number_bins, modal_intensity + N * index_bin))
        background_noise = bin_edges[thresh]
    elif N==0:
        background_noise = bin_edges[modal_intensity]

    SNR = np.max(wave_env) - background_noise
    SN = np.array([frame-background_noise-activity_threshold_dB for frame in wave_env])
    acoustic_activity = np.sum([i > 0 for i in SN])/float(len(SN))


    # Compute acoustic events
    start_event = [n[0] for n in np.argwhere((SN[:-1] < 0) & (SN[1:] > 0))]
    end_event = [n[0] for n in np.argwhere((SN[:-1] > 0) & (SN[1:] < 0))]
    if len(start_event)!=0 and len(end_event)!=0:
        if start_event[0]<end_event[0]:
            events=list(zip(start_event, end_event))
        else:
            events=list(zip(end_event, start_event))
        count_acoustic_events = len(events)
        average_duration_e = np.mean([end - begin for begin,end in events] )
        average_duration_s = average_duration_e * file.duration / float(len(SN))
    else:
        count_acoustic_events = 0
        average_duration_s = 0


    dict = {'SNR' : SNR, 'Acoustic_activity' : acoustic_activity, 'Count_acoustic_events' : count_acoustic_events, 'Average_duration' : average_duration_s}
    return dict




#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def remove_noiseInSpectro(spectro, histo_relative_size=8, window_smoothing=5, N=0.1, dB=False, plot=False):
    """

    Compute a new spectrogram which is "Noise Removed".

    spectro: spectrogram of the audio signal
    histo_relative_size: ration between the size of the spectrogram and the size of the histogram
    window_smoothing: number of points to apply a mean filtering on the histogram and on the background noise curve
    N: Parameter to set the threshold around the modal intensity
    dB: If set at True, the spectrogram is converted in decibels
    plot: if set at True, the function plot the orginal and noise removed spectrograms

    Output:
        Noise removed spectrogram

    Ref: Towsey, Michael W. (2013) Noise removal from wave-forms and spectrograms derived from natural recordings of the environment.
    Towsey, Michael (2013), Noise Removal from Waveforms and Spectrograms Derived from Natural Recordings of the Environment. Queensland University of Technology, Brisbane.
    """

    low_value = 1.e-07 # Minimum value for the new spectrogram (preferably slightly higher than 0)
    half_window_smoothing = int(window_smoothing /2)

    if dB:
        spectro = 20*np.log10(spectro)

    len_spectro_e = len(spectro[0])
    histo_size = int(len_spectro_e/histo_relative_size)

    background_noise=[]
    for row in spectro:
        hist, bin_edges = np.histogram(row, bins=histo_size, density=False)

        hist_smooth = ([np.mean(hist[i - half_window_smoothing: i + half_window_smoothing]) for i in range(half_window_smoothing, len(hist) - half_window_smoothing)])
        hist_smooth = np.concatenate((np.zeros(half_window_smoothing), hist_smooth, np.zeros(half_window_smoothing)))


        modal_intensity = int(np.min([np.argmax(hist_smooth), 95 * histo_size / 100])) # test if modal intensity value is in the top 5%

        if N>0:
            count_thresh = 68 * sum(hist_smooth) / 100
            count = hist_smooth[modal_intensity]
            index_bin = 1
            while count < count_thresh:
                if modal_intensity + index_bin <= len(hist_smooth):
                    count = count + hist_smooth[modal_intensity + index_bin]
                if modal_intensity - index_bin >= 0:
                    count = count + hist_smooth[modal_intensity - index_bin]
                index_bin += 1
            thresh = int(np.min((histo_size, modal_intensity + N * index_bin)))
            background_noise.append(bin_edges[thresh])
        elif N==0:
            background_noise.append(bin_edges[modal_intensity])


    background_noise_smooth = ([np.mean(background_noise[i - half_window_smoothing: i + half_window_smoothing]) for i in range(half_window_smoothing, len(background_noise) - half_window_smoothing)])
    # keep background noise at the end to avoid last row problem (last bin with old microphones)
    background_noise_smooth = np.concatenate((background_noise[0:half_window_smoothing], background_noise_smooth, background_noise[-half_window_smoothing:]))

    new_spec = np.array([col - background_noise_smooth for col in spectro.T]).T
    new_spec = new_spec.clip(min=low_value) # replace negative values by value close to zero

    #Figure
    if plot:
        colormap="jet"
        fig = plt.figure()
        a = fig.add_subplot(1,2,1)
        if dB:
            plt.imshow(new_spec, origin="lower", aspect="auto", cmap=colormap, interpolation="none")
        else:
            plt.imshow(20*np.log10(new_spec), origin="lower", aspect="auto", cmap=colormap, interpolation="none")
        a = fig.add_subplot(1,2,2)
        if dB:
            plt.imshow(new_spec, origin="lower", aspect="auto", cmap=colormap, interpolation="none")
        else:
            plt.imshow(20*np.log10(spectro), origin="lower", aspect="auto", cmap=colormap, interpolation="none")
        plt.show()



    return new_spec


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_NB_peaks(spectro, frequencies, freqband = 200, normalization= True, slopes=(0.01,0.01)):
    """

    Counts the number of major frequency peaks obtained on a mean spectrum.

    spectro: spectrogram of the audio signal
    frequencies: list of the frequencies of the spectrogram
    freqband: frequency threshold parameter (in Hz). If the frequency difference of two successive peaks is less than this threshold, then the peak of highest amplitude will be kept only.
    normalization: if set at True, the mean spectrum is scaled between 0 and 1
    slopes: amplitude slope parameter, a tuple of length 2. Refers to the amplitude slopes of the peak. The first value is the left slope and the second value is the right slope. Only peaks with higher slopes than threshold values will be kept.

    Ref: Gasc, A., Sueur, J., Pavoine, S., Pellens, R., & Grandcolas, P. (2013). Biodiversity sampling using a global acoustic approach: contrasting sites with microendemics in New Caledonia. PloS one, 8(5), e65311.

    """

    meanspec = np.array([np.mean(row) for row in spectro])

    if normalization:
         meanspec =  meanspec/np.max(meanspec)

    # Find peaks (with slopes)
    peaks_indices = np.r_[False, meanspec[1:] > np.array([x + slopes[0] for x in meanspec[:-1]])] & np.r_[meanspec[:-1] > np.array([y + slopes[1] for y in meanspec[1:]]), False]
    peaks_indices = peaks_indices.nonzero()[0].tolist()

    #peaks_indices = signal.argrelextrema(np.array(meanspec), np.greater)[0].tolist() # scipy method (without slope)


    # Remove peaks with difference of frequency < freqband
    nb_bin=next(i for i,v in enumerate(frequencies) if v > freqband) # number of consecutive index
    for consecutiveIndices in [np.arange(i, i+nb_bin) for i in peaks_indices]:
        if len(np.intersect1d(consecutiveIndices,peaks_indices))>1:
            # close values has been found
            maxi = np.intersect1d(consecutiveIndices,peaks_indices)[np.argmax([meanspec[f] for f in np.intersect1d(consecutiveIndices,peaks_indices)])]
            peaks_indices = [x for x in peaks_indices if x not in consecutiveIndices] # remove all inddices that are in consecutiveIndices
            peaks_indices.append(maxi) # append the max
    peaks_indices.sort()


    peak_freqs = [frequencies[p] for p in peaks_indices] # Frequencies of the peaks

    return len(peaks_indices)

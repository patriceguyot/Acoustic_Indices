#!/usr/bin/env python

"""
    Set of functions to compute acoustic indices in the framework of Soundscape Ecology.

    Some features are inspired or ported from those proposed in:
        - seewave R package (http://rug.mnhn.fr/seewave/) / Jerome Sueur, Thierry Aubin and  Caroline Simonis
        - soundecology R package (http://cran.r-project.org/web/packages/soundecology/index.html) / Luis J. Villanueva-Rivera and Bryan C. Pijanowski

    This file use an object oriented type for audio files described in the file "acoustic_index.py".



"""

__author__ = "Patrice Guyot"
__version__ = "0.2"
__credits__ = ["Patrice Guyot", "Alice Eldridge", "Mika Peck"]
__email__ = ["guyot.patrice@gmail.com", "alicee@sussex.ac.uk", "m.r.peck@sussex.ac.uk"]
__status__ = "Development"


from scipy import signal
import numpy as np



#--------------------------------------------------------------------------------------------------------------------------------------------------
def compute_spectrogram(file, windowLength=512, windowHop= 256, scale_audio=True, square=True, windowType='hanning', centered=False, normalized = False ):
    """
    Compute a spectrogram of an audio signal.
    Return a list of list of values as the spectrogram, and a list of frequencies.

    Keyword arguments:
    file -- the real part (default 0.0)

    Parameters:
    @param sig: the audio signal samples (read for example from a .wav with scipy.io.wavfile)
    wLen: length of the fft window (in samples)
    wHop: hop size of the fft window (in samples)
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

    if centered:
        time_shift = int(windowLength/2)
        times = range(time_shift, len(sig)+1-time_shift, windowHop) # centered
        frames = [sig[i-time_shift:i+time_shift]*W for i in times] # centered frames
    else:
        times = range(0, len(sig)-windowLength+1, windowHop)
        frames = [sig[i:i+windowLength]*W for i in times]

    if square:
        spectro =  [abs(np.fft.rfft(frame, windowLength))[0:windowLength/2]**2 for frame in frames]
    else:
        spectro =  [abs(np.fft.rfft(frame, windowLength))[0:windowLength/2] for frame in frames]

    spectro=np.transpose(spectro) # set the spectro in a friendly way

    if normalized:
        spectro = spectro/np.max(spectro) # set the maximum value to 1 y

    frequencies = [e * file.niquist / float(windowLength / 2) for e in range(windowLength / 2)] # vector of frequency<-bin in the spectrogram
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
    In theis code, the Bioacoustic Index correspond to the area under the mean spectre (in dB) minus the minimum frequency value of this mean spectre.

    Reference: Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne remote sensing. Ecological Applications 17: 2137-2144.

    spectro: the spectrogram of the audio signal
    min_freq: minimum frequency (in Hertz)
    max_freq: maximum frequency (in Hertz)

    Ported from the soundecology R package.
    """

    min_freq_bin=np.argmin([abs(e - min_freq) for e in frequencies]) # min freq in samples (or bin)
    max_freq_bin=np.ceil(np.argmin([abs(e - max_freq) for e in frequencies])) # max freq in samples (or bin)

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

    env = abs(signal.hilbert(sig)) # Modulo of the Hilbert Envelope

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

    frequencies, pxx = signal.welch(file.sig_float, fs=file.sr, window='hamming', nperseg=windowLength, noverlap=windowLength/2, nfft=windowLength, detrend=False, return_onesided=True, scaling='density', axis=-1) # Estimate power spectral density using Welch's method
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

    spectro: spectrogram of the audio signal
    freq_band_Hz: frequency band (in Hertz)
    max_freq: the maximum frequency to consider (in Hertz)
    db_threshold: the minimum dB value to consider for the bins of the spectrogram
    freq_step: width of the frequency bands that are considered

    Ported from the soundecology R package.
    """


    #TODO better comments for frequencies bands



    bands_Hz = range(0, max_freq, freq_step)
    bands_bin = [f / freq_band_Hz for f in bands_Hz]

    spec_AEI = 20*np.log10(spectro/np.max(spectro))
    spec_AEI_bands = [spec_AEI[bands_bin[k]:bands_bin[k]+bands_bin[1],] for k in range(len(bands_bin))]

    values = [np.sum(spec_AEI_bands[k]>db_threshold)/float(spec_AEI_bands[k].size) for k in range(len(bands_bin))]

    return gini(values)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_ADI(spectro, freq_band_Hz,  max_freq=10000, db_threshold=-50, freq_step=1000):
    """
    Compute Acoustic Diversity Index.

    spectro: spectrogram of the audio signal
    freq_band_Hz: frequency band (in Hertz)
    max_freq: the maximum frequency to consider (in Hertz)
    db_threshold: the minimum dB value to consider for the bins of the spectrogram
    freq_step: width of the frequency bands that are considered


    Ported from the soundecology R package.
    """

     #TODO better comments for frequencies bands

    bands_Hz = range(0, max_freq, freq_step)
    bands_bin = [f / freq_band_Hz for f in bands_Hz]

    spec_AEI = 20*np.log10(spectro/np.max(spectro))
    spec_AEI_bands = [spec_AEI[bands_bin[k]:bands_bin[k]+bands_bin[1],] for k in range(len(bands_bin))]

    values = [np.sum(spec_AEI_bands[k]>db_threshold)/float(spec_AEI_bands[k].size) for k in range(len(bands_bin))]

    # Shannon Entropy of the values
    #shannon = - sum([y * np.log(y) for y in values]) / len(values)  # Follows the R code. But log is generally log2 for Shannon entropy. Equivalent to shannon = False in soundecology.

    # The following is equivalent to shannon = True (default) in soundecology. Compute the Shannon diversity index from the R function diversity {vegan}.
    #v = [x/np.sum(values) for x in values]
    #v2 = [-i * j  for i,j in zip(v, np.log(v))]
    #return np.sum(v2)

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
    Compute the spectral centroid of an audio signal from its spectrogram
    """

    centroid = [ np.sum(magnitudes*frequencies) / np.sum(magnitudes) for magnitudes in spectro.T]


    return centroid

#!/usr/bin/env python

"""
    compute_indices.py
    Set of functions to compute acoustic indices in the framework of Soundscape Ecology.

    Some features are inspired or ported from those proposed in:
        - seewave R package (http://rug.mnhn.fr/seewave/) / Jerome Sueur, Thierry Aubin and  Caroline Simonis
        - soundecology R package (https://cran.r-project.org/web/packages/soundecology/index.html) / Luis J. Villanueva-Rivera and Bryan C. Pijanowski

"""

__author__ = "Patrice Guyot"
__version__ = "0.1"
__credits__ = ["Patrice Guyot", "Alice Eldridge", "Mika Peck"]
__email__ = ["guyot.patrice@gmail.com", "alicee@sussex.ac.uk", "m.r.peck@sussex.ac.uk"]
__status__ = "Development"



from scipy import signal
import numpy as np



def pcm2float(sig, dtype='float64'):
    """Convert PCM signal to floating point with a range from -1 to 1.

    Use dtype='float32' for single precision.

    Parameters
    ----------
    sig : array_like
        Input array, must have integral type.
    dtype : data type, optional
        Desired (floating point) data type.

    Returns
    -------
    numpy.ndarray
        Normalized floating point data.

    See Also
    --------
    float2pcm, dtype

    """
    sig = np.asarray(sig)
    if sig.dtype.kind not in 'iu':
        raise TypeError("'sig' must be an array of integers")
    dtype = np.dtype(dtype)
    if dtype.kind != 'f':
        raise TypeError("'dtype' must be a floating point type")

    i = np.iinfo(sig.dtype)
    abs_max = 2 ** (i.bits - 1)
    offset = i.min + abs_max
    return (sig.astype(dtype) - offset) / abs_max


#--------------------------------------------------------------------------------------------------------------------------------------------------
def compute_spectrogram(sig, wLen=512, wHop= 256, scale_audio=True, square=True, windowType='hanning', centered=False, normalization = False ):
    """
    Compute spectrogram of an audio signal.
    Return a list of list of values.

    sig: the audio signal samples (read for example from a .wav with scipy.io.wavfile)
    wLen: length of the fft window (in samples)
    wHop: hop size of the fft window (in samples)
    scale_audio: if set as True, the signal samples are scale between -1 and 1 (as the audio convention). If false the signal samples remains Integers (as output from scipy.io.wavfile)
    square: if set as True, the spectrogram is computed as the square of the magnitude of the fft. If not, it is the magnitude of the fft.
    hamming: if set as True, the spectrogram use a correlation with a hamming window.
    centered: if set as true, each resulting fft is centered on the corresponding sliding window
    """

    if scale_audio:
        sig = pcm2float(sig, 'float64') # transform signal int to float between -1 and 1


    W = signal.get_window(windowType, wLen, fftbins=False)

    if centered:
        time_shift = int(wLen/2)
        times = range(time_shift, len(sig)+1-time_shift, wHop) # centered
        frames = [sig[i-time_shift:i+time_shift]*W for i in times] # centered frames
    else:
        times = range(0, len(sig)-wLen+1, wHop)
        frames = [sig[i:i+wLen]*W for i in times]

    if square:
        spectro =  [abs(np.fft.rfft(frame, wLen))[0:wLen/2]**2 for frame in frames]
    else:
        spectro =  [abs(np.fft.rfft(frame, wLen))[0:wLen/2] for frame in frames]

    spectro=np.transpose(spectro) # set the spectro in a friendly way

    if normalization:
        spectro = spectro/np.max(spectro)

    return spectro





#-----------------------------------------------------------------------------
def compute_ACI(spectro,j_bin):
    """
    Compute the Acoustic Complexity Index from the spectrogram of an audio signal.
    spectro: the spectrogram of the audio signal
    j_bin: temporal size of the frame (in samples)

    Ported from soundecology R package.
    """

    #times = range(0, spectro.shape[1], j_bin) # relevant time indices
    times = range(0, spectro.shape[1]-10, j_bin) # alternative time indices to follow the R code

    jspecs = [np.array(spectro[:,i:i+j_bin]) for i in times]  # sub-spectros of temporal size j

    aci = [sum((np.sum(abs(np.diff(jspec)), axis=1) / np.sum(jspec, axis=1))) for jspec in jspecs] 	# list of ACI values on each jspecs
    return sum(aci)





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_BI(spectro, frequencies, min_freq = 2000, max_freq = 8000):
    """
    Compute the Bioacoustic Index from the spectrogram of an audio signal.
    spectro: the spectrogram of the audio signal
    min_freq: minimum frequency (in Hertz)
    max_freq: maximum frequency (in Hertz)

    Ported from the soundecology R package.
    """


    min_freq_bin=np.argmin([abs(e - min_freq) for e in frequencies]) # min freq in samples (or bin)
    max_freq_bin=np.ceil(np.argmin([abs(e - max_freq) for e in frequencies])) # max freq in samples (or bin)

    min_freq_bin = min_freq_bin - 1 # alternative value to follow the R code



    spec_BI = 20 * np.log10(spectro/np.max(spectro))  #  Use of decibel values. Equivalent in the R code to: spec_left <- spectro(left, f = samplingrate, wl = fft_w, plot = FALSE, dB = "max0")$amp
    spec_BI_mean = 10 * np.log10 (np.mean(10 ** (spec_BI/10), axis=1))     # Equivalent in the R code to: return(a*log10(mean(10^(x/a))))
    spec_BI_mean_segment =  spec_BI_mean[min_freq_bin:max_freq_bin]
    spec_BI_mean_segment_normalized = spec_BI_mean_segment - min(spec_BI_mean_segment)
    area = np.sum(spec_BI_mean_segment_normalized / (frequencies[1]-frequencies[0]))   # Equivalent in the R code to: left_area <- sum(specA_left_segment_normalized * rows_width)

    return (area)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_SH(spectro):
    """
    Compute Spectral Entropy of Shannon from the spectrogram of an audio signal.
    spectro: the spectrogram of the audio signal

    Ported from the seewave R package.
    """
    N = spectro.shape[0]
    spec = np.sum(spectro,axis=1)
    spec = spec / np.sum(spec)  # Normalization
    # TODO check the normalization (why not in the spectogram)
    return - sum([y * np.log2(y) for y in spec]) / np.log2(N)  #Equivalent in the R code to: z <- -sum(spec*log(spec))/log(N)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_TH(sig):
    """
    Compute Temporal Entropy of Shannon from an audio signal.

    Ported from the seewave R package.
    """
    env = abs(signal.hilbert(sig)) # Modulo of the Hilbert Envelope

    env = env / np.sum(env)  # Normalization
    N = len(env)
    return - sum([y * np.log2(y) for y in env]) / np.log2(N)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_NDSI(sig, sr, winLen = 1024, anthrophony=[1000,2000], biophony=[2000,11000]):
    """
    Compute Normalized Difference Sound Index from an audio signal.
    anthrophony: list of two values containing the minimum and maximum frequencies (in Hertz) for antrophony
    biophony: list of two values containing the minimum and maximum frequencies (in Hertz) for biophony

    Inspired by the seewave R package, the soundecology R package and the original matlab code from the authors.
    """

    #TODO put winLen in the call of the function

    frequencies, pxx = signal.welch(sig, fs=sr, window='hamming', nperseg=winLen, noverlap=winLen/2, nfft=winLen, detrend=False, return_onesided=True, scaling='density', axis=-1) # Estimate power spectral density using Welch's method
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

    Inspired by http://mathworld.wolfram.com/GiniCoefficient.html and https://en.wikipedia.org/wiki/Gini_coefficient
    """
    y = sorted(values)
    n = len(y)
    G = np.sum([i*j for i,j in zip(y,range(1,n+1))])
    G = 2 * G / np.sum(y) - (n+1)
    return G/n


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_AEI(sig, sr,  max_freq=10000, db_threshold=-50, freq_step=1000):
    """
    Compute Acoustic Evenness Index of an audio signal.
    sig: the audio signal
    sr: the sampling rate of the audio signal
    max_freq: the maximum frequency to consider (in Hertz)
    db_threshold: the minimum dB value to consider for the bins of the spectrogram
    freq_step: width of the frequency bands that are considered

    Ported from the soundecology R package.
    """

    #TODO better comments for frequencies bands
    freq_band_Hz = max_freq / freq_step
    windowLength = sr / freq_band_Hz
    spectro = compute_spectrogram(sig, wLen=windowLength, wHop= windowLength, scale_audio=True, square=False, windowType='hanning', centered=False, normalization= False )
    #TODO why not considering the spectro outside of the function ?

    bands_Hz = range(0, max_freq, freq_step)
    bands_bin = [f / freq_band_Hz for f in bands_Hz]

    spec_AEI = 20*np.log10(spectro/np.max(spectro))
    spec_AEI_bands = [spec_AEI[bands_bin[k]:bands_bin[k]+bands_bin[1],] for k in range(len(bands_bin))]

    values = [np.sum(spec_AEI_bands[k]>db_threshold)/float(spec_AEI_bands[k].size) for k in range(len(bands_bin))]

    return gini(values)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_ADI(sig, sr,  max_freq=10000, db_threshold=-50, freq_step=1000):
    """
    Compute Acoustic Diversity Index.
    sig: the audio signal
    sr: the sampling rate of the audio signal
    max_freq: the maximum frequency to consider (in Hertz)
    db_threshold: the minimum dB value to consider for the bins of the spectrogram
    freq_step: width of the frequency bands that are considered


    Ported from the soundecology R package.
    """

    freq_band_Hz = max_freq / freq_step
    windowLength = sr / freq_band_Hz
    spectro = compute_spectrogram(sig, wLen=windowLength, wHop= windowLength, scale_audio=True, square=False, windowType='hanning', centered=False, normalization= False )
    #TODO why not considering the spectro outside of the function ?

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
def compute_zcr(sig, wLen=512, wHop= 256):
    """
    Compute the Zero Crossing Rate of an Audio signal
    """
    times = range(0, len(sig)-wLen+1, wHop)
    frames = [sig[i:i+wLen] for i in times]
    return [len(np.where(np.diff(np.signbit(x)))[0])/float(wLen) for x in frames]


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def compute_rms_energy(sig, wLen=512, wHop=256):
    """
    Compute the RMS short time energy
    """
    times = range(0, len(sig)-wLen+1, wHop)
    frames = [sig[i:i+wLen] for i in times]
    return [np.sqrt(sum([ x**2 for x in frame ]) / wLen) for frame in frames]



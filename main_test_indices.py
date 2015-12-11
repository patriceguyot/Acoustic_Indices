__author__ = 'guyot'

#!/usr/bin/env python

"""
    Compute and output acoustic indices
"""

__author__ = "Patrice Guyot"
__version__ = "0.1"
__credits__ = ["Patrice Guyot", "Alice Eldridge", "Mika Peck"]
__email__ = ["guyot.patrice@gmail.com", "alicee@sussex.ac.uk", "m.r.peck@sussex.ac.uk"]
__status__ = "Development"



from scipy.io.wavfile import read as wavread
from compute_indices import *
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------------------
def plot_spectrogram(data,sig,sr):
    """
    Plot spectrogram
    """
    colormap="jet"
    npdata=np.array(data)
    ims = 20.*np.log10(np.abs(npdata)/10e-6) # amplitude to decibel
    plt.imshow(ims, origin="lower", aspect="auto", cmap=colormap, interpolation="none")
    plt.xlabel("time ")
    plt.ylabel("frequency")

    binsize=2**10
    timebins, freqbins = np.shape(ims)
    xlocs = np.float32(np.linspace(0, timebins-1, 5))
    plt.xticks(xlocs, ["%.02f" % l for l in ((xlocs*len(sig)/timebins)+(0.5*binsize))/sr])


    plt.show()




if __name__ == '__main__':


    # Read signal -------------------------------------


    filename = 'PL-16_0_20150605_0615.wav'
    filename = '1.wav'


    print 'Read signal ', filename
    sr, sig = wavread('data/' + filename)

    f_nyquist = sr/2


    # -----------------------------------------------------
    print('Compute spectrogram')
    windowLength = 512  # Length of sliding window for fft (samples)
    windowHop = 512     # Length of lag for fft (samples)

    frequencies = [e * f_nyquist / float(windowLength / 2) for e in range(windowLength / 2)] # vector of frequency<-bin in the spectrogram
    spectro = compute_spectrogram(sig, wLen=windowLength, wHop= windowHop, scale_audio=True, square=False, windowType='hanning', centered=False, normalization= False )
    plot_spectrogram(spectro, sig, sr)

    # -----------------------------------------------------
    print('RMS energy')
    windowLength = 512  # Length of teh frame (samples)
    windowHop = 256     # Length of the lag (samples)

    rms = compute_rms_energy(pcm2float(sig, 'float64'), windowLength, windowHop)
    plt.plot(rms)
    plt.show()

    # -----------------------------------------------------
    print('ZCR')
    windowLength = 512  # Length of teh frame (samples)
    windowHop = 256     # Length of the lag (samples)

    zcr = compute_zcr(sig, windowLength, windowHop)
    plt.plot(zcr)
    plt.show()

    # -----------------------------------------------------
    print('Acoustic Diversity Index')
    adi =  compute_ADI(sig, sr, max_freq = 10000, db_threshold = -50, freq_step = 1000)
    print(adi)

    # -----------------------------------------------------
    print('Acoustic Evenness Index')
    aei =  compute_AEI(sig, sr, max_freq = 10000, db_threshold = -50, freq_step = 1000)
    print(aei)


    # -----------------------------------------------------
    print('NDSI')
    windowLength = 1024
    #spectro = compute_spetrogram(sig, wLen=windowLength, wHop= windowHop, scale_audio=True, square=True, windowType='hanning', centered=False, normalization= True )
    # TODO  check the pcm2float signal int to float between -1 and 1
    ndsi =  compute_NDSI(pcm2float(sig, 'float64'), sr, windowLength, [1000,2000], [2000,11000])
    print(ndsi)

    # -----------------------------------------------------
    print('Spectral Entropy')
    spectro = compute_spectrogram(sig, wLen=windowLength, wHop= windowHop, scale_audio=True, square=False, windowType='hanning', centered=False, normalization= False )
    sh = compute_SH(spectro)
    print(sh)

    # -----------------------------------------------------
    print('Compute Temporal Entropy')
    th = compute_TH(sig)
    print(th)

    # -----------------------------------------------------
    print('Compute Total Entropy')
    H = sh * th
    print(H)

    # -----------------------------------------------------
    print('Compute ACI')
    spectro = compute_spectrogram(sig, wLen=windowLength, wHop= windowHop, scale_audio=False, square=False, windowType='hamming', centered=False, normalization= True )

    j = 5 # j in seconds
    j_bin = j * sr / windowHop
    aci = compute_ACI(spectro, j_bin)
    print(aci)

    # -----------------------------------------------------
    print('Compute BI')
    spectro = compute_spectrogram(sig, wLen=windowLength, wHop= windowHop, scale_audio=True, square=False, windowType='hanning', centered=False, normalization= False )
    min_freq = 2000
    max_freq = 8000
    BI = compute_BI(spectro, frequencies, min_freq, max_freq)
    print(BI)

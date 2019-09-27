__author__ = 'guyot'

#!/usr/bin/env python

"""
    Compute and output acoustic indices
"""

__author__ = "Patrice Guyot"
__version__ = "0.4"
__credits__ = ["Patrice Guyot", "Alice Eldridge", "Mika Peck"]
__email__ = ["guyot.patrice@gmail.com", "alicee@sussex.ac.uk", "m.r.peck@sussex.ac.uk"]
__status__ = "Development"


from compute_indice import *
from acoustic_index import *
import yaml
from scipy import signal
import csv



if __name__ == '__main__':

    #Set config file
    yml_file = 'yaml/config_014_butter.yaml'

    # Read signal -------------------------------------
    filename = 'audio_files/BALMER-02_0_20150620_0445.wav'


    file = AudioFile(filename, verbose=True)

    with open(yml_file, 'r') as stream:
        data_config = yaml.load(stream, Loader=yaml.FullLoader)


    # Pre-processing -----------------------------------------------------------------------------------
    if 'Filtering' in data_config:
        if data_config['Filtering']['type'] == 'butterworth':
            print('- Pre-processing - High-Pass Filtering:', data_config['Filtering'])
            freq_filter = data_config['Filtering']['frequency']
            Wn = freq_filter/float(file.niquist)
            order = data_config['Filtering']['order']
            [b,a] = signal.butter(order, Wn, btype='highpass')
            # to plot the frequency response
            #w, h = signal.freqz(b, a, worN=2000)
            #plt.plot((file.sr * 0.5 / np.pi) * w, abs(h))
            #plt.show()
            file.process_filtering(signal.filtfilt(b, a, file.sig_float))
        elif data_config['Filtering']['type'] == 'windowed_sinc':
            print('- Pre-processing - High-Pass Filtering:', data_config['Filtering'])
            freq_filter = data_config['Filtering']['frequency']
            fc = freq_filter / float(file.sr)
            roll_off = data_config['Filtering']['roll_off']
            b = roll_off / float(file.sr)
            N = int(np.ceil((4 / b)))
            if not N % 2: N += 1  # Make sure that N is odd.
            n = np.arange(N)
            # Compute a low-pass filter.
            h = np.sinc(2 * fc * (n - (N - 1) / 2.))
            w = np.blackman(N)
            h = h * w
            h = h / np.sum(h)
            # Create a high-pass filter from the low-pass filter through spectral inversion.
            h = -h
            h[(N - 1) / 2] += 1
            file.process_filtering(np.convolve(file.sig_float, h))



    # Compute Indices -----------------------------------------------------------------------------------
    print('- Compute Indices')
    ci = data_config['Indices'] # use to simplify the notation
    for index_name in ci:  # iterate over the index names (key of dictionary in the yml file)


        if index_name == 'Acoustic_Complexity_Index':
            print('\tCompute', index_name)
            spectro, _ = compute_spectrogram(file, **ci[index_name]['spectro'])
            methodToCall = globals().get(ci[index_name]['function'])
            j_bin = int(ci[index_name]['arguments']['j_bin'] * file.sr / ci[index_name]['spectro']['windowHop']) # transform j_bin in samples
            main_value, temporal_values = methodToCall(spectro, j_bin)
            file.indices[index_name] = Index(index_name, temporal_values=temporal_values, main_value=main_value)


        elif index_name == 'Acoustic_Diversity_Index':
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            freq_band_Hz = ci[index_name]['arguments']['max_freq'] / ci[index_name]['arguments']['freq_step']
            windowLength = int(file.sr / freq_band_Hz)
            spectro,_ = compute_spectrogram(file, windowLength=windowLength, windowHop= windowLength, scale_audio=True, square=False, windowType='hanning', centered=False, normalized= False )
            main_value = methodToCall(spectro, freq_band_Hz, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Acoustic_Evenness_Index':
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            freq_band_Hz = ci[index_name]['arguments']['max_freq'] / ci[index_name]['arguments']['freq_step']
            windowLength = int(file.sr / freq_band_Hz)
            spectro,_ = compute_spectrogram(file, windowLength=windowLength, windowHop= windowLength, scale_audio=True, square=False, windowType='hanning', centered=False, normalized= False )
            main_value = methodToCall(spectro, freq_band_Hz, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Bio_acoustic_Index':
            print('\tCompute', index_name)
            spectro, frequencies = compute_spectrogram(file, **ci[index_name]['spectro'])
            methodToCall = globals().get(ci[index_name]['function'])
            main_value = methodToCall(spectro, frequencies, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Normalized_Difference_Sound_Index':
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            main_value = methodToCall(file, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'RMS_energy':
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            temporal_values = methodToCall(file, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, temporal_values=temporal_values)


        elif index_name == 'Spectral_centroid':
            print('\tCompute', index_name)
            spectro, frequencies = compute_spectrogram(file, **ci[index_name]['spectro'])
            methodToCall = globals().get(ci[index_name]['function'])
            temporal_values = methodToCall(spectro, frequencies)
            file.indices[index_name] = Index(index_name, temporal_values=temporal_values)


        elif index_name == 'Spectral_Entropy':
            print('\tCompute', index_name)
            spectro, _ = compute_spectrogram(file, **ci[index_name]['spectro'])
            methodToCall = globals().get(ci[index_name]['function'])
            main_value = methodToCall(spectro)
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Temporal_Entropy':
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            main_value = methodToCall(file, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'ZCR':
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            temporal_values = methodToCall(file, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, temporal_values=temporal_values)


        elif index_name == 'Wave_SNR':
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            values = methodToCall(file, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, values=values)


        elif index_name == 'NB_peaks':
            print('\tCompute', index_name)
            spectro, frequencies = compute_spectrogram(file, **ci[index_name]['spectro'])
            methodToCall = globals().get(ci[index_name]['function'])
            main_value = methodToCall(spectro, frequencies, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Acoustic_Diversity_Index_NR': # Acoustic_Diversity_Index with Noise Removed spectrograms
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            freq_band_Hz = ci[index_name]['arguments']['max_freq'] / ci[index_name]['arguments']['freq_step']
            windowLength = int(file.sr / freq_band_Hz)
            spectro,_ = compute_spectrogram(file, windowLength=windowLength, windowHop= windowLength, scale_audio=True, square=False, windowType='hanning', centered=False, normalized= False )
            spectro_noise_removed = remove_noiseInSpectro(spectro, **ci[index_name]['remove_noiseInSpectro'])
            main_value = methodToCall(spectro_noise_removed, freq_band_Hz, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Acoustic_Evenness_Index_NR': # Acoustic_Evenness_Index with Noise Removed spectrograms
            print('\tCompute', index_name)
            methodToCall = globals().get(ci[index_name]['function'])
            freq_band_Hz = ci[index_name]['arguments']['max_freq'] / ci[index_name]['arguments']['freq_step']
            windowLength = int(file.sr / freq_band_Hz)
            spectro,_ = compute_spectrogram(file, windowLength=windowLength, windowHop= windowLength, scale_audio=True, square=False, windowType='hanning', centered=False, normalized= False )
            spectro_noise_removed = remove_noiseInSpectro(spectro, **ci[index_name]['remove_noiseInSpectro'])
            main_value = methodToCall(spectro_noise_removed, freq_band_Hz, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Bio_acoustic_Index_NR': # Bio_acoustic_Index with Noise Removed spectrograms
            print('\tCompute', index_name)
            spectro, frequencies = compute_spectrogram(file, **ci[index_name]['spectro'])
            spectro_noise_removed = remove_noiseInSpectro(spectro, **ci[index_name]['remove_noiseInSpectro'])
            methodToCall = globals().get(ci[index_name]['function'])
            main_value = methodToCall(spectro_noise_removed, frequencies, **ci[index_name]['arguments'])
            file.indices[index_name] = Index(index_name, main_value=main_value)


        elif index_name == 'Spectral_Entropy_NR': # Spectral_Entropy with Noise Removed spectrograms
            print('\tCompute', index_name)
            spectro, _ = compute_spectrogram(file, **ci[index_name]['spectro'])
            spectro_noise_removed = remove_noiseInSpectro(spectro, **ci[index_name]['remove_noiseInSpectro'])
            methodToCall = globals().get(ci[index_name]['function'])
            main_value = methodToCall(spectro_noise_removed)
            file.indices[index_name] = Index(index_name, main_value=main_value)




    # Output Indices -----------------------------------------------------------------------------------
    print('- Write Indices')
    writer = csv.writer(open('dict.csv', 'w'))

    keys = ['filename']
    values = [file.file_name]
    for index, Index in file.indices.items():
        for key, value in Index.__dict__.items():
            if key != 'name':
                keys.append(index + '__' + key)
                values.append(value)
    writer.writerow(keys)
    writer.writerow(values)

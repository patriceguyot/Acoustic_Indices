#!/usr/bin/python

__author__ = 'guyot'
__version__ = "0.4"


import numpy as np
from scipy.io.wavfile import read as wavread
from scipy.io.wavfile import write as wavwrite
from os.path import basename


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

#------------------------------------------------------------------------
def float2pcm(sig, dtype='int16'):
    """Convert floating point signal with a range from -1 to 1 to PCM.

    Any signal values outside the interval [-1.0, 1.0) are clipped.
    No dithering is used.

    Note that there are different possibilities for scaling floating
    point numbers to PCM numbers, this function implements just one of
    them.  For an overview of alternatives see
    http://blog.bjornroche.com/2009/12/int-float-int-its-jungle-out-there.html

    Parameters
    ----------
    sig : array_like
        Input array, must have floating point type.
    dtype : data type, optional
        Desired (integer) data type.

    Returns
    -------
    numpy.ndarray
        Integer data, scaled and clipped to the range of the given
        `dtype`.

    See Also
    --------
    pcm2float, dtype

    """
    sig = np.asarray(sig)
    if sig.dtype.kind != 'f':
        raise TypeError("'sig' must be a float array")
    dtype = np.dtype(dtype)
    if dtype.kind not in 'iu':
        raise TypeError("'dtype' must be an integer type")

    i = np.iinfo(dtype)
    abs_max = 2 ** (i.bits - 1)
    offset = i.min + abs_max
    return (sig * abs_max + offset).clip(i.min, i.max).astype(dtype)



#------------------------------------------------------------------------
class AudioFile(object):

    def __init__(self, file_path, verbose=False):

        default_channel = 0  # Used channel if the audio file contains more than one channel
        #default_channel = 1

        if verbose:
            print('Read the audio file:', file_path)
        try:
            sr, sig = wavread(file_path)
        except IOError:
            print("Error: can\'t read the audio file:", file_path)
        else:
            if verbose:
                print('\tSuccessful read of the audio file:', file_path)
            if len(sig.shape)>1:
                if verbose:
                    print('\tThe audio file contains more than one channel. Only the channel', default_channel, 'will be used.')
                sig=sig[:, default_channel] # Assign default channel
            self.sr = sr
            self.sig_int = sig
            self.sig_float = pcm2float(sig,dtype='float64')
            self.niquist = sr/2
            self.file_path = file_path
            self.file_name = basename(file_path)
            self.filtered = False
            self.duration = len(sig)/float(sr)
            self.indices = dict()  # empty dictionary of Index

    def process_filtering(self, sig_float_filtered, write=False, output_file_name = None):
        self.filtered = True
        self.sig_int = float2pcm(sig_float_filtered)
        self.sig_float = sig_float_filtered
        if write:
            wavwrite(output_file_name, file.sr, file.sig_int)

#----------------------------------------------------------------------------------------
class Index(object):

    def __init__(self, name, temporal_values=None, main_value=None, values=None):
        self.name = name
        if main_value is not None:
            self.main_value = main_value
        if temporal_values is not None:
            self.min = np.nanmin(temporal_values)
            self.max = np.nanmax(temporal_values)
            self.mean = np.nanmean(temporal_values)
            self.median = np.nanmedian(temporal_values)
            self.std = np.nanstd(temporal_values)
            self.var = np.nanvar(temporal_values)
        if values is not None:
            for name,value in values.items():
                setattr(self, name, value)

    def print_all(self):
        print("Name:", self.name)
        if hasattr(self, 'main_value'):
            print("Main value:", self.main_value)
        if hasattr(self, 'temporal_values'):
            print("Minimum:", self.min)
            print("Maximum:", self.max)
            print("Arithmetic mean value:", self.mean)
            print("Median value:", self.median)
            print("Standard deviation:", self.std)
            print("Variance:", self.var)
        if hasattr(self, 'values'):
            for name,value in self.values.iteritems():
                print(name, ':', value)

if __name__ == '__main__':
    print("Hello")
    centroid = Index('centroid',2,[2,3,4,5,19,10])
    centroid.print_all()

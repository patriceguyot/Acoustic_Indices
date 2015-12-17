# Acoustic_Indices

Acoustic_Indices is a Python library to extract acoustic indices from a set of audio files, in the framework of Soundscape Ecology.


## Prerequisites

 * [Numpy](http://www.numpy.org/)
 * [Scipy](http://www.scipy.org/)
 * [Matlplotlib](http://matplotlib.org/) (for graphing)
 * [PyYAML](http://pyyaml.org/wiki/PyYAMLDocumentation) (to read configuration file)

## Features and indices

* Read WAV files (using scipy)
* Features extraction from Soundscape Ecology 
    * Acoustic Complexity Index
    * Acoustic Diversity Index
    * Acoustic Evenness Index
    * Bioacoustic Index
    * Normalized Difference Sound Index
    * Spectral Entropy
    * Temporal Entropy
    
* Spectral features extraction
    * Spectral centroid 
    * Spectrogram
    
* Temporal features extraction
    * RMS energy
    * Zero Crossing Rate
    

## Usage

$python compute_indices.py

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

Versions: 

* 0.2: yaml configuration file. Object oriented audio file and index.
* 0.1: First commit 



## Credits

TODO: Write credits

## License

TODO: Write license


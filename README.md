# Acoustic_Indices

Acoustic_Indices is a Python library to extract global acoustic indices from an audio file for use as a biodiversity proxy, within the framework of Ecoacoustics.


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
    * Number of Peaks
    * Wave Signal to Noise Ratio

* Spectral features extraction
    * Spectral centroid
    * Spectrogram
    * Noise removed spectrogram

* Temporal features extraction
    * RMS energy
    * Zero Crossing Rate


## Usage

$python main\_test\_indices.py

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

Versions:

* 0.4: Port from Python 2 to Python 3
* 0.3: New features: wave SNR, spectro noise removed, NB_peaks.
* 0.2: yaml configuration file. Object oriented audio file and index.
* 0.1: First commit



## Credits

The following indices are based on the following papers and inspired in part by the R packages [seewave] (https://cran.r-project.org/package=seewave) and [soundecology] (https://cran.r-project.org/package=soundecology)  
* Acoustic Complexity Index - Pieretti et al. (2011)
* Acoustic Diversity Index - Villanueva-Rivera et al. (2011)
* Acoustic Evenness Index - Villanueva-Rivera et al. (2011)
* Bioacoustic Index - Boelman, et al. (2007)
* Normalized Difference Sound Index - Kasten et al. (2012)
* Spectral Entropy - Sueur et al. (2008)
* Temporal Entropy - Sueur et al. (2008)

Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne remote sensing. Ecological Applications 17: 2137-2144.

Farina A, Pieretti N, Piccioli L (2011) The soundscape methodology for long-term bird monitoring: a Mediterranean Europe case-study. Ecological Informatics, 6, 354-363.

Kasten, E.P., Gage, S.H., Fox, J. & Joo, W. (2012). The remote environmental assessment laboratory's acoustic library: an archive for studying soundscape ecology. Ecological Informatics, 12, 50-67.

Pieretti N, Farina A, Morri FD (2011) A new methodology to infer the singing activity of an avian community: the Acoustic Complexity Index (ACI). Ecological Indicators, 11, 868-873.

Sueur, J., Pavoine, S., Hamerlynck, O. & Duvail, S. (2008) - Rapid acoustic survey for biodiversity appraisal. PLoS ONE, 3(12): e4065.

Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246. doi: 10.1007/s10980-011-9636-9.


This research was generously funded by Leverhulme Research Project Grant RPG-2014-403

## License

GPL V 3 - see [LICENSE.txt](https://github.com/patriceguyot/Acoustic_Indices/blob/master/LICENSE.txt) for more info

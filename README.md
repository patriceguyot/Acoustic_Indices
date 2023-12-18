# Acoustic_Indices

Acoustic_Indices is a Python library to extract global acoustic indices from an audio file for use as a biodiversity proxy, within the framework of Ecoacoustics.


## Indices

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




## Reference

If you use this code, please cite: Alice Eldridge, & Patrice Guyot. (2023). Python implementation of acoustic indices and low level descriptors (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.10391651


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10391651.svg)](https://doi.org/10.5281/zenodo.10391651)





## Publications

This code have been used in the following scientific papers:

* Durbridge, S., & Murphy, D. T. (2023). Assessment of soundscapes using self-report and physiological measures. Acta Acustica, 7, 6.
* Martínez-Tabares, F., & Orozco-Alzate, M. (2022, November). Selection of acoustic features for the discrimination between highly and moderately transformed Colombian soundscapes. In Ibero-American Conference on Artificial Intelligence (pp. 121-132). Cham: Springer International Publishing.
* Sumitani, S., Suzuki, R., Morimatsu, T., Matsubayashi, S., Arita, T., Nakadai, K., & Okuno, H. G. (2020, January). Soundscape Analysis of Bird Songs in Forests Using Microphone Arrays. In 2020 IEEE/SICE International Symposium on System Integration (SII) (pp. 634-639). IEEE.
* Carruthers-Jones, J., Eldridge, A., Guyot, P., Hassall, C., & Holmes, G. (2019). The call of the wild: Investigating the potential for ecoacoustic methods in mapping wilderness areas. Science of the Total Environment, 695, 133797.
* Eldridge, A., Guyot, P., Moscoso, P., Johnston, A., Eyre-Walker, Y., & Peck, M. (2018). Sounding out ecoacoustic metrics: Avian species richness is predicted by acoustic indices in temperate but not tropical habitats. Ecological Indicators, 95, 939-952.



## Usage

Test that everything is going well on one audio file:

``` 
$python main\_test\_indices.py 
```

Compute indices from a directory of audio files:
```
$python  main\_compute\_indices\_from\_dir
```

This is the pipeline of the processing:

* Get list of indices and features from a Yaml configuration file 
* Read WAV files (using scipy)
* For each audio file, compute stats (min, max, mean, median, std, var) for temporal acoustic indices or get global value for other indices
* Output a csv file (with a row for each audio file and a column for each index)

## Prerequisites

 * [Numpy](http://www.numpy.org/)
 * [Scipy](http://www.scipy.org/)
 * [Matlplotlib](http://matplotlib.org/) (for graphing)
 * [PyYAML](http://pyyaml.org/wiki/PyYAMLDocumentation) (to read configuration file)





## History

Versions:

* 0.5: New main (main\_compute\_indices\_from\_dir) to compute all indices from a directory of audio files
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


This research was generously funded by Leverhulme Research Project Grant RPG-2014-403.





# A biological reading of a palimpsest
Laura C. Viñas-Caron*, Ismael Rodríguez Palomo*, Natasha Fazlic, Jiří Vnouček, Matthew Driscoll, Sarah Fiddyment, Matthew J. Collins

Code for the MS data analysis for the paper "A biological reading of a palimpsest".

MALDI-TOF data is deposited in Zenodo with DOI: https://doi.org/10.5281/zenodo.6967158

The code is divided in two notebooks:
- Spectra preprocessing from sample replicate spectra to a feature matrix of samples by peaks. We do square root intensity transformation, smoothing, baseline correction, and total ion current (TIC) normalization. We then detect peaks using the SuperSmoother for noise estimation. Peaks are aligned by calculating warping functions using LOWESS and binned to adjust their mass numerically and uncommon peaks within replicates removed. Finally, replicates are merged by averaging their peak intensities and the feature matrix is generated.

- Data analysis. Batch correction with surrogate variable analysis, PCA, clustering, binary discriminant analysis and theoretical collagen markers alignment using bacollite

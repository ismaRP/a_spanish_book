# A biological reading of a palimpsest
Laura C. Viñas-Caron*, Ismael Rodríguez Palomo*, Natasha Fazlic, Jiří Vnouček, Matthew Driscoll, Sarah Fiddyment, Matthew J. Collins

Code for the MS data analysis for the paper "A biological reading of a palimpsest".

The code is divided in two notebooks:
- Spectra preprocessing. From sample replicate spectra to a matrix of samples by peaks. We do square root intensity transformation, smoothing, baseline correction, and total ion current (TIC) normalization. We then detect peaks using the SuperSmoother for noise estimation. Peaks are aligned by calculating warping functions using LOWESS. Finally peaks are binned to finally adjust them numerically and uncommon peaks within replicates removed.

- Data analysis. Batch correction with surrogate variable analysis, PCA, clustering, binary discriminant analysis and theoretical collagen markers alignment using bacollite

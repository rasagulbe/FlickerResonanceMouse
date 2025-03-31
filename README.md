# FlickerResonanceMouse

## The code accompanying the paper "Spatiotemporal resonance in mouse primary visual cortex"

### Rasa Gulbinaite, Mojtaba Nazari, Michael E. Rule, Edgar J. Bermudez-Contreras, Mike X Cohen, Majid H. Mohajerani, J. Alexander Heimel, 2024, _Current Biology_, 34, 4184-4196.e7
https://www.sciencedirect.com/science/article/abs/pii/S0960982224010595
##
The repository contains MATLAB code for:
1. Alignment of recordings across days and to Allen Common Coordinate Framework (CCFv3). 
2. Computation of trial-average and single-trial power and phase maps.
3. Preprocessing raw widefield imaging data (loading, alignment between recordings, global signal regression (GSR), baseline correction, robust regression, spatial Gaussian filter).
4. Testing statistical significance of responses to flicker in V1 and whole cortex level (bootstrapping using Wilcox's bootstrap-t).
5. Spatiotemporal pattern analysis using two-stage GED: (1) dimensionality reduction using PCA; (2) generalized eigenvalue decomposition (GED) as implemented in Rhythmic Source Separation method (RESS; see Cohen and Gulbinaite, 2017 NeuroImage; PMID: 27916666).

#### NOTE: All raw data can be found on DRYAD repository (https://datadryad.org/dataset/doi:10.5061/dryad.vdncjsz42)

**Requirements:** MATLAB Image Processing Toolbox and Signal Processing Toolbox. The code was run on MATLAB 2019b.

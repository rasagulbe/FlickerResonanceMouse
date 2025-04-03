### The code accompanying the publication
## [Spatiotemporal resonance in mouse primary visual cortex](https://www.cell.com/current-biology/abstract/S0960-9822(24)01059-5)

### Rasa Gulbinaite, Mojtaba Nazari, Michael E. Rule, Edgar J. Bermudez-Contreras, Mike X Cohen, Majid H. Mohajerani, J. Alexander Heimel, 2024, _Current Biology_, 34, 4184-4196.e7
#### DOI: [https://www.cell.com/current-biology/abstract/S0960-9822(24)01059-5](https://www.cell.com/current-biology/abstract/S0960-9822(24)01059-5)
##
![Graphical abstract Git](https://github.com/user-attachments/assets/2351831d-18cc-4eb8-99e0-1daa31e28eea)
##
### The repository contains MATLAB code for:

1. Alignment of recordings across days and to Allen Common Coordinate Framework (CCFv3). 
2. Computation of trial-average and single-trial power and phase maps.
3. Preprocessing raw widefield imaging data (loading, alignment between recordings, global signal regression (GSR), baseline correction, robust regression, spatial Gaussian filter).
4. Testing statistical significance of responses to flicker in V1 and whole cortex level (bootstrapping using Wilcox's bootstrap-t).
5. Spatiotemporal pattern analysis using two-stage GED:
     * dimensionality reduction using PCA;
     * generalized eigenvalue decomposition (GED) as implemented in Rhythmic Source Separation method (RESS; see Cohen and Gulbinaite, 2017 NeuroImage; [DOI: https://doi.org/10.1016/j.neuroimage.2016.11.036 ](https://doi.org/10.1016/j.neuroimage.2016.11.036)).
##
### Dataset 
All raw data can be found on DRYAD repository (https://datadryad.org/dataset/doi:10.5061/dryad.vdncjsz42)
##
### Requirements
MATLAB Image Processing Toolbox and Signal Processing Toolbox. The code was run on MATLAB 2019b.

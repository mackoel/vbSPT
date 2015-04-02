MLD: diffusion constant estimators for 1-state diffusive SPT data with 
localization errors and blur.

ML1_logLlambda.m - log likelihood function for camera-based single
particle tracking data with camera-based blur and localization errors
that are given for every point.

This code makes use of the routine triSym_d1Inv_trjWise.mex (source
code in HMMcore) to compute the determinant and some elements of the
inverse of a symmetric tridiagonal matrix.

vestergaardCVE.m - Covariance-based estimator of diffusion constant
and (average) localization errors, by Vestergaard, Blainey, and
Flyvbjerg. Consistency checks and confidence interval not (yet)
implemented.

lsqMSD.m - an msd-based estimator (note that the CVE estimator
generally has better performance).

ML1_preprocess_mixed_columns : data preprocessor.

Martin Lindén, bmelinden@gmail.com, 2015-04-02

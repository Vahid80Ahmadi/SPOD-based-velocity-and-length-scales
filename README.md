# Spectral characteristics and a new Strouhal number scaling of noise emitting open hydrogen-enriched swirl-stabilized turbulent premixed flames

The spectraly resolved velocity and length scales of noise sources in an open hydrogen-enriched swril-stablized can be defined using Spectral Proper Orthagonal Decomposition (SPOD). 

This code includes the required MATALB\textsuperscript{\textregistered} scripts for calculating the velocity and lenght scales, and is cosidered as a part of supplementary material of [1]. The mathematical and physical explanation of this code is provided in [1]. 

- The *velocity_scale.m* and *lenght_scale.m* use the most energetic SPOD mode at each frequency to calculate velocity and lenght scales, respectively. 

- To obtain the SPOD modes, *main_spod.m* is used. The SPOD calculation was done using *spod.m*. For more details about SPOD please refer to [2].

- The *spod.m* function can be either found in this [link](https://www.mathworks.com/matlabcentral/fileexchange/65683-spectral-proper-orthogonal-decomposition-spod) or this [link](https://github.com/SpectralPOD/spod_matlab) as well.

- The *license.txt* summarizes the terms and conditions for *spod.m*.

- For local wavenumber calculations, a 2D regularized differentiation method similar to [3] was used. The relevant functions are included in *velocity_scale.m*. Since the included functions are just a MATLAB\textsuperscript{\textregistered} version of the python codes, which can be found in this [link](https://github.com/rickchartrand/regularized_differentiation/tree/master), provided by the author of [3], *lisence2.txt* explains the terms and conditions of the mentioned original code.

- Raw data will be provided upon request. 

- If you are using this code to calculate length and velocity scales, please kindly cite [1].

- If you have any question, please do not hesitate to contact us via either of these emails.
  - Sina.kheirkhah@ubc.ca
  - vahid.ahmadi@ubc.ca
  - ahmadi.v.1380@gmail.com

# References

[1] V. Ahmadi, J. Fleger, J. Ho, S. Mohammadnejad, M. Talei, S. Kheirkhah, Spectral characteristics and a new Strouhal number scaling of noise emitting open hydrogen-enriched swirl-stabilized turbulent premixed flames, Int. J. Heat Fluid Flow.	

[2]  A. Towne, O. T. Schmidt, T. Colonius, Spectral proper orthogonal decomposition and its relationship to dynamic mode decomposition and resolvent analysis, J. of Fluid Mech. 847, 821–867, 2018.

[3] R. Chartrand, Numerical differentiation of noisy, nonsmooth, multidimensional data, in: 2017 IEEE Global Conference on Signal and Information Processing (GlobalSIP), IEEE, 2017, pp. 244–248.

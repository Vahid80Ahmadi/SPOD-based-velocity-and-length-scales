# Spectral characteristics and a new Strouhal number scaling of noise emitting open hydrogen-enriched swirl-stabilized turbulent premixed flames

The spectraly resolved velocity and length scales of noise sources in an open hydrogen-enriched swril-stablized can be defined using Spectral Proper Orthagonal Decomposition (PSOD). 

This code includes the required MATALB scripts for calculating the velocity and lenght scales, and is cosidered as a part of supplementary material of [1]. The mathematical and physical explanation of this code is provided in [1]. 

- The *velocity_scale.m* and *lenght_scale.m* use the most energetic SPOD mode at each frequency to calculate velocity and lenght scales, respectively. 

- To obtain the SPOD modes, *main_spod.m* is used. The SPOD calculation was done using *spod.m*. For more details about SPOD please refer to [2].

- The *spod.m* function can be find in this [link](https://www.mathworks.com/matlabcentral/fileexchange/65683-spectral-proper-orthogonal-decomposition-spod) or this [link](https://github.com/SpectralPOD/spod_matlab) as well.

- The *license.txt* summarize the terms and conditions for *spod.m*.

- If you are using this code to calculate length and velocity scales, please kindly cite [1].

- If you have any question, please do not hesitate to contact us via either of these emails.
  - Sina.kheirkhah@ubc.ca
  - vahid.ahmadi@ubc.ca
  - ahmadi.v.1380@gmail.com

Refrences:

***[1]*** V. Ahmadi, J. Fleger, J. Ho, S. Mohammadnejad, M. Talei, S. Kheirkhah, Spectral characteristics and a new Strouhal number scaling of noise emitting open hydrogen-enriched swirl-stabilized turbulent premixed flames, Int. J. Heat Fluid Flow.	

***[2]*** Towne, A., Schmidt, O. T., Colonius, T., Spectral proper orthogonal decomposition and its relationship to dynamic mode decomposition and resolvent analysis, J. of Fluid Mech. 847, 821–867, 2018

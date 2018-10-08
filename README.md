## HFR CS PROCESSING TOOLBOX ##

Copyright (C) 2018 Brian Emery, Ph.D


Tools for processing oceanographic HF radar cross spectra with direction
finding methods. 


WHAT THIS DOES
- Provides research quality software for processing HFR data following the
  methods of Lipe et al. (2006).
- Employs a home made ship removal algorithm that follows what is known 
  about how it's really done. 
- Generalizes the cross spectra data structure for arbitrary arrays.
- Allows the use of the imageFOL toolbox by Anthony Kirincich 
  (https://github.com/akirincich/imageFOLs).
- Includes several DF methods used in Emery (2018), and an estimate of 
  MUSIC error from Stoica and Nehorai (1989) used in Emery and Washburn (2018).
- Uses my version of CODAR's method for single vs dual determination. This 
  part is a work in progress for other arrays. 


HOW TO USE IT
- start with run_cs_processing.m 


TO DO
- I need to add a demo/test, e.g. using Anthony's data for validation, but
  also could build some edge case tests using simiulations (eg 0-360 
  transition, etc)


ACKNOWLEDGMENT
- AK, HFRP, COS


REFERENCES
Lipa, B., B. Nyden, D. S. Ullman, and E. Terrill, (2006). SeaSonde radial 
  velocities: Derivation and internal consistency.
  IEEE Journal of Oceanic Engineering, 31 (4) 850?861.
Stoica, P., & Nehorai, A. (1989). MUSIC, maximum likelihood, and 
  Cramer-Rao bound. IEEE Transactions on Acoustics, Speech, and Signal 
  Processing, 37(5), 720-741.
Emery, B. and Washburn, L. (2018). Uncertainty Estimates for SeaSonde HF 
  Radar Ocean Current Observations. Journal of Atmospheric and Oceanic 
  Technology, Submitted.
Emery, B. (2018). Evaluation of Alternative Direction of Arrival Methods
  for Oceanographic HF Radars. IEEE Journal of Oceanic Engineering, 
  Submitted.


ARCHITECTURE GOALS
- could follow range processing
- arbitrary array geometry, fft length, etc
- arbitrary doa method 
- future use of things like different detectors, MAP algorithm, ...
- tests using Anthony's data for validation, but also could build some
  edge case tests using simiulations (eg 0-360 transition, etc)

NOTES
- .m files use functions as blocks of code - code folding (cmd =) makes it easy to move among these.
- Data structures contain variables of similar origin following the HFRProgs
  convention, rows = locations, cols = time
- Data structures are initialized with appropriately named function 
  (e.g. doa_struct.m) to enable standardization


CODING PRINCIPLES
(The goal, if not the reality)
- Minimize repetition (dont repeat yourself)
- Make code re-usable and recyclable. Make general functions. 
- Code should be pretty and readable sentences and paragraphs, aid the reader when
  possible
- balance between future usages (flexibility) and getting the current job done (purpose built)
- Good design is simple
- Write computer programs to make them easy for people to read.
- Have  a clear division between code that is custom for a particular application, 
  and the general/easily repurposed code





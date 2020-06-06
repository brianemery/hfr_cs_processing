## HFR CS PROCESSING TOOLBOX FOR MATLAB ##

v1.0

[![DOI](https://zenodo.org/badge/84593561.svg)](https://zenodo.org/badge/latestdoi/84593561)

Tools for processing oceanographic HF radar cross spectra with direction
finding methods. 

- Provides research quality software for processing HFR data following the
  methods of Lipa et al. (2006).
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
- download it and add the unzipped directory to your MATLAB path
- run the demo to make sure it works (run_cs_processing_demo.m)
- edit run_cs_processing.m for more advanced applications 
- edit doa_on_cs.m (~line 80) to include the use of parfor 



TO DO
- Include installation instructions?
- Could build some edge case tests using simiulations (eg 0-360 
  transition, etc)
- Arbitrary arrays need a detection method (that is, single bearing vs dual
  vs ... etc).
- This commit includes many more mfiles than are actually needed due to MATLAB's 
  dependency tool - need a better way to isolate tools. 


ACKNOWLEDGMENT

The release is self contained but includes code from the following people
and or toolboxes. HFRProgs [1] by David Kaplan set the standard, and data
structures here follow the same basic formatting. I've include a few mfiles
from HFRprogs as dependencies. As mentioned above, this toolbox can use 
imageFOL by Anthony Kirincich [2]. Finally, the release includes code obtained
from CODAR Ocean Sensors for reading cross spectra files. 

[1] https://github.com/rowg/hfrprogs  
[2] https://github.com/akirincich/imageFOLs   


REFERENCES

Lipa, B., B. Nyden, D. S. Ullman, and E. Terrill, (2006). SeaSonde radial   
&nbsp;&nbsp;velocities: Derivation and internal consistency.  
&nbsp;&nbsp;IEEE Journal of Oceanic Engineering, 31 (4) 850?861.  

Stoica, P., & Nehorai, A. (1989). MUSIC, maximum likelihood, and  
&nbsp;&nbsp;Cramer-Rao bound. IEEE Transactions on Acoustics, Speech, and Signal     
&nbsp;&nbsp;Processing, 37(5), 720-741.  

Emery, B. and Washburn, L. (2019). Uncertainty Estimates for SeaSonde HF   
&nbsp;&nbsp;Radar Ocean Current Observations. Journal of Atmospheric and 
&nbsp;&nbsp;Oceanic Technology 36.2: 231-247.

Emery, B. (2019). Evaluation of Alternative Direction of Arrival Methods  
&nbsp;&nbsp;for Oceanographic HF Radars. IEEE Journal of Oceanic Engineering,   
&nbsp;&nbsp;doi: 10.1109/JOE.2019.2914537.  


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

- Minimize repetition (dont repeat yourself)
- Make code re-usable and recyclable. Make general functions. 
- Code should be pretty and readable sentences and paragraphs, aid the reader when
  possible
- balance between future usages (flexibility) and getting the current job done (purpose built)
- Good design is simple
- Write computer programs to make them easy for people to read.
- Have  a clear division between code that is custom for a particular application, 
  and the general/easily repurposed code
  
  
HOW TO CITE WITH BIBTEX
  
```
@misc{Emery2018code,
  author       = {Brian Emery},
  title        = {{HFR CS Processing Toolbox for MATLAB}, Software Release Version 1.0, https://doi.org/10.5281/zenodo.1451950},
  version      = {1.0},
  month        = oct,
  year         = 2018,
  doi          = {10.5281/zenodo.1451950},
  url          = {https://doi.org/10.5281/zenodo.1451950}
}
```

VERSION NOTES

1.5 (5 June 2020)
Updates include a demonstration (run_cs_processing_demo.m) and test data from BML1. This update also includes the likelihood ratio 
detection method from an in-prep manuscript (suitable for use with MLE direction finding). 



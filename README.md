## HFR CS PROCESSING TOOLBOX FOR MATLAB ##

v2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5598294.svg)](https://doi.org/10.5281/zenodo.5598294)

Tools for processing oceanographic HF radar cross spectra with direction
finding methods. 

- Provides research quality software for processing HFR data following the
  methods of Lipa et al. (2006).
- Employs a home made ship removal algorithm that follows what is known 
  about how it's really done. (Presently disabled) 
- Generalizes the cross spectra data structure for arbitrary arrays.
- Allows the use of the imageFOL toolbox by Anthony Kirincich 
  (https://github.com/akirincich/imageFOLs).
- Includes several DF methods used in Emery (2018), and an estimate of 
  MUSIC error from Stoica and Nehorai (1989) used in Emery and Washburn (2018).
- Uses my version of CODAR's method for single vs dual determination. Can use 
  the Likelihood Ratio (See Emery et. al 2022) for larger arrays. 
- Will work with data from LERA and WERA. It has been used to 
  process LERA data.


HOW TO USE IT
- download it and cd to the unzipped directory
- run install_hfr_cs_proc.m to modify the MATLAB path
- run the demo to make sure it works (run_cs_processing_demo.m)
- edit run_cs_processing.m for more advanced applications 
- edit doa_on_cs.m (~line 80) to include the use of parfor 


TO DO
- Could build some edge case tests using simiulations (eg 0-360 
  transition, etc)
- This commit includes many more mfiles than are actually needed due to MATLAB's 
  dependency tool - need a better way to isolate tools. 


ACKNOWLEDGMENT

The release is self contained but includes code from the following people
and or toolboxes. HFRProgs [1] by David Kaplan set the standard, and data
structures here follow the same basic formatting. I've include a few mfiles
from HFRprogs as dependencies. As mentioned above, this toolbox can use 
imageFOL by Anthony Kirincich [2]. The toolbox includes code obtained
from CODAR Ocean Sensors for reading cross spectra files. Version 1.5 
includes data from BML1 provided by William Speiser and John Largier. 

[1] https://github.com/rowg/hfrprogs  
[2] https://github.com/akirincich/imageFOLs   


REFERENCES

Lipa, B., B. Nyden, D. S. Ullman, and E. Terrill, (2006). SeaSonde radial velocities: Derivation and internal consistency. IEEE Journal of Oceanic Engineering, 31 (4) 850?861.  

Stoica, P., & Nehorai, A. (1989). MUSIC, maximum likelihood, and the Cramer-Rao bound. IEEE Transactions on Acoustics, Speech, and Signal Processing, 37(5), 720-741.  

Emery, B. and Washburn, L. (2019). Uncertainty Estimates for SeaSonde HF Radar Ocean Current Observations. Journal of Atmospheric and Oceanic Technology 36.2: 231-247.

Emery, B. (2019). Evaluation of Alternative Direction of Arrival Methods for Oceanographic HF Radars. IEEE Journal of Oceanic Engineering, doi: 10.1109/JOE.2019.2914537.  

Emery, B., Kirincich, A., and Washburn, L. (2022). Direction Finding and Likelihood Ratio Detection for 
Oceanographic HF Radars. Journal of Atmospheric and Oceanic Technology, 39(2), 223-235. https://journals.ametsoc.org/view/journals/atot/39/2/JTECH-D-21-0110.1.xml  

ARCHITECTURE GOALS

- arbitrary array geometry, fft length, etc
- arbitrary doa method 


NOTES

- .m files use functions as blocks of code - code folding (cmd =) makes it easy to move among these.
- Data structures contain variables of similar origin following the HFRProgs
  convention, rows = locations, cols = time
- Data structures are initialized with appropriately named function 
  (e.g. doa_struct.m) to enable standardization
- End product is presently a structure that I'm calling a DOA structure, which is roughly
  equivalent to a SeaSonde short-time radial (pre-merge), plus a lot more meta data.


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
  
  
HOW TO CITE

(See the Zenodo site to cite specific older versions): 

Emery, B., 2021: HFR CS Processing Toolbox for MATLAB, software release version 2.0, https://doi.org/10.5281/zenodo.5598294. URL: https://doi.org/10.5281/zenodo.5598294, doi:10.4445281/zenodo.5598294.

Using Bibtex:
 
```
@misc{Emery2021code,
  author       = {Brian Emery},
  title        = {{HFR CS Processing Toolbox for MATLAB}, Software Release Version 2.0, https://doi.org/10.5281/zenodo.5598294},
  version      = {2.0},
  month        = oct,
  year         = 2021,
  doi          = {10.5281/zenodo.5598294},
  url          = {https://doi.org/10.5281/zenodo.5598294}
}
```

VERSION NOTES

1.0 (18 July 2018)
Original version can be found here: [![DOI](https://zenodo.org/badge/84593561.svg)](https://zenodo.org/badge/latestdoi/84593561)

1.5 (16 June 2020)
Updates include an install script (install_hfr_cs_proc.m) and demonstration code (run_cs_processing_demo.m) that 
uses test data from BML1. This update also includes the likelihood ratio detection method from an in-prep
manuscript (suitable for use with MLE direction finding - more about this at a later time). 

2.0 (21 Oct 2021) 
Updates, improvements and new features related to Emery et al. 2022.

2.1 (1 June 2022)
Updates including radial metrics, many code fixes and likely bug introductions. 

## HFR CS PROCESSING TOOLBOX ##

Copyright (C) 2017 Brian Emery


ARCHITECTURE IDEAS
- could follow range processing
- arbitrary array geometry, fft length, etc
- arbitrary doa method 
- future use of things like different detectors, MAP algorithm, ...
- tests using Anthony's data for validation, but also could build some
  edge case tests using simiulations (eg 0-360 transition, etc)

  
NOTES

Main calling function is run_cs_processing.m

Dependencies (most of them) are in /private directory. These may include HFRProgs and 
other common stuff

.m files use functions as blocks of code - code folding (cmd =) makes it easy to move 
among these.


DATA STRUCTURES

- Contain variables of similar scope, rows = locations, cols = time
- Initialized with appropriately named function (e.g. doa_struct.m) to enable 
  standardization


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





function APM = apm_struct(tf)
% APM_STRUCT Initialize a standard APM structure
% APM = apm_struct
%
% Initialize the APM output structure.
%
% Use these:
% BEAR A13R A13I A23R A23I FLAG A13M A13P A23M A23P 
%
% for field names (they come out of loop files)
%
% Specify/standardize that fields are column vectors
%
% Given the input 'true', creates an expanded APM struct for APMs derived
% from ships of opportunity.
%
% Also when the array matrix is specified in the APM struct, orient it 
% #bearings x #elements, while the output of get_array_matrix is output as
% #elements x #bearings, eg:
%     APM.A = APM.A.';
%
% NOTE ABOUT UNITS
% APMs seem to usually be in 'degCWN', but RADIALstruct has
% Degrees_ccw_from_east. In run_cs_processing.m, get_radial_meta makes the
% change from deg cwN to ccwE. 
%
% SEE ALSO
% apm_struct_ext.m 

% Copyright (C) 2009-2010 Brian M. Emery
% June 2009

% TO DO
% add this to extended:
% targe distance in meters
% S.TRGD
% might want to use other LOOP file keys as field names, see
% write_loop_files.m

% optionally expand struct for AIS derived patterns
if nargin < 1, tf = 0; end

% General, HFRP like stuff
APM.Type =  'Measured Antenna Pattern' ;
APM.Creator = 'B. Emery';
APM.SiteName   = [];
APM.SiteOrigin = [];
APM.FileName   = [];

% Time
APM.CreateTimeStamp = datestr(now);


% define standard fields
vars = {'README','BEAR','A13R','A13RQ','A13I','A13IQ', ...
                        'A23R','A23RQ','A23I','A23IQ', ...
                        'A13M','A13P' ,'A23M' ,'A23P', ...
                        'A33R','A33I','Units','A'};
                   
% create empties
for i = 1:numel(vars)
        APM.(vars{i}) = [];               
end

APM.README.A = 'Array Matrix, #bearings x #elements';

%------------------------------
%  EXPANSION FOR AIS derived APM 
%--------------------------------
if tf
   % expand README
   APM.README.BEAR          = 'Bearing from HF site to ship using AIS';
   APM.README.BEARMean      = 'Mean of AIS for each CSQ time';
   APM.README.BEARCOS       = 'Bearing of ship using MUSIC and measured APM [cwN]';
   APM.README.fbinIdx       = 'Freq Bin Index in the CSQ';
   APM.README.dopplerRadVel = 'Doppler spectrum bin velocities of ship peak';
   APM.README.meanAPM       = 'Mean of APM values derived from each CSQ ship peak';


   vars = {'SNR','BEARMean','BEARCOS','fbinIdx','dopplerRadVel', ...
       'rCellIdx'};
   
   % create empties
   for i = 1:numel(vars)
       APM.(vars{i}) = [];
   end

   % Empty cell for csq names
   APM.csqFileName={};

end
    

APM.ProcessingSteps = {};


end

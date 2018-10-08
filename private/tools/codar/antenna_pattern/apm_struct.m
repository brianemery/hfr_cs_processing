function APM = apm_struct(tf)
%APM_STRUCT Initialize a standard APM structure
% Initialize the APM output structure.
%
% Use these:
% BEAR A13R A13I A23R A23I FLAG A13M A13P A23M A23P 
%
% for field names (they come out of loop files)
%
% Specify/standardize that fields are column vectors

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
                        'A33R','A33I','Units'};
                   
% create empties
for i = 1:numel(vars)
        APM.(vars{i}) = [];               
end


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
    
% 
% % Obsolete version:
% return
% [APM.README, ...
% APM.a13r, ...
% APM.a13i, ...
% APM.a23r, ...
% APM.a23i, ...
% APM.a33r, ...
% APM.a33i, ...
% APM.a13m, ...
% APM.a13p, ...
% APM.a23m, ...
% APM.a23p, ...
% APM.snr, ...
% APM.bear, ...
% APM.bearMean, ... 
% APM.bearCOS, ...
% APM.fbinIdx, ...
% APM.dopplerRadVel, ...
% APM.rCellIdx] = deal([]);
% 
% APM.csqFileName={};
% 

APM.ProcessingSteps = {};


end

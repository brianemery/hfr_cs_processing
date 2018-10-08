function APM = apm_struct_ext(N)
% APM STRUCT EXT - extended APM struct initialization
% APM = apm_struct_ext(N)
%
% Initialize the APM output structure. *Pre-allocates and
% includes non-standard fields* used in the calculation of 
% APMs from ships of opportunity.


% Copyright (C) 2009-2010 Brian M. Emery
% June 2009


% General, HFRP like stuff
APM.Type = ['Measured Antenna Pattern'];
APM.Creator = 'B. Emery';
APM.SiteName   = 'XXXX';
APM.SiteOrigin = [NaN NaN];
APM.FileName   = [];


%  EXPANSION FOR AIS derived APM
% expand README
APM.README.TimeStamp     = 'Interpolated from AIS based on RadVel';
APM.README.BEAR          = 'Bearing from HF site to ship using AIS';
APM.README.fbinIdx       = 'Freq Bin Index in the CSQ';
APM.README.dopplerRadVel = 'Doppler spectrum bin velocities of ship peak (cm/s)';
APM.README.SNR           = 'Various SNR calculations (snr_calculations.m)';
APM.README.shipID        = 'Ship ID from AIS';
APM.README.Evalue        = 'Eigen Values sm to lg';
APM.README.Evec1         = 'Noise Eigen Vector';
APM.README.Evec2         = 'Noise Eigen Vector';
APM.README.Evec3         = 'Signal Eigen Vector (APM components)';
APM.README.RadVel        = 'Mean Ship Radial Velocity and stats'; 
APM.README.Lon           = 'Interpolated from AIS based on CS radvel';
APM.README.Lat           = 'Interpolated from AIS based on CS radvel';


% DEFINE STANDARD FIELDS
vars = {'TimeStamp','BEAR', 'Lon','Lat'...
        'A13R','A13I','A23R','A23I', ...
        'A13M','A13P','A23M' ,'A23P', ...
        'A33R','A33I'};
APM = create_empties(APM,vars,N,1);


% DEFINE NON-STANDARD FIELDS
vars = {'dopplerRadVel','shipID'...
        'fbinIdx','rCellIdx'};
APM = create_empties(APM,vars,N,1);


% MEAN SHIP RADIAL VELOCITY AND STATS
APM.RadVel.xbar =[];
vars={'xbar','stdev','hi','lo','median','N'};
APM.RadVel = create_empties(APM.RadVel,vars,N,1);


% EIGEN VALUES AND VECTORS
vars = {'Evalue','Evec1','Evec2','Evec3'};
APM = create_empties(APM,vars,N,3);

% SNR DATA
APM.SNR(1).README = {'stdCodar  =  Standard Codar Method'; ...
                     'backGrnd  =  Above Backgound (from cs_filter.m)' ; ...
                     'localDopp =  Local adjacent doppler bins (no bragg)' ; ...
                     'inRange   =  Adjacent range cells'};
vars = {'stdCodar','backGrnd','localDopp','inRange'};
APM.SNR = create_empties(APM.SNR,vars,N,3);


% Empty cells
APM.csqFileName = cell(N,1);
APM.ProcessingSteps = {};

% Time, other meta info
APM.CreationInfo = creation_info;


end
% --------------------------
function APM = create_empties(APM,vars,r,c)

for i = 1:numel(vars)
   APM.(vars{i}) = NaN(r,c);
end


end

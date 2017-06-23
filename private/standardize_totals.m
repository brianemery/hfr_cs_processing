function TUV = standardize_totals(AngDiff,U,V,gridd,serialtime,Err)
% STANDARDIZE TOTALS - Put tot_*.mat file variables in a TUV structure
% 
% TUV = standardize_totals(AngDiff,U,V,gridd,serialtime,Err)
%
% Probably will have lots of use for this in converting from my older
% format to HFRP format for total data

% Copyright (C) 2009-2010 Brian Emery
% 19 Nov 2009

if nargin < 6, 
    disp('standardize_totals.m: need to include Err'), Err = [];
end

% create empty TUV stucture:
TUV = TUVstruct([1 1],1);

% Fill in some information
TUV.DomainName = 'UCSB';
TUV.CreationInfo = 'BME';
TUV.CreateTimeStamp = now;

% Fill in data
TUV.LonLat = gridd;
TUV.TimeStamp = serialtime;
TUV.U = U;
TUV.V = V;
TUV.AngDiff = AngDiff;

% add details
TUV.OtherMetadata.AngDiff = 'Angle closest to 90deg in least square fit to radials';
TUV.ProcessingSteps{end+1} = mfilename;

% Fit difference
% Note that I think these are different than DMK's rms fit diff
TUV.ErrorEstimates.Type = 'Mean Fit Difference (see totcalc.m)';
TUV.ErrorEstimates.Err = Err;
TUV.ErrorEstimates.TotalErrors = Err; % For HFRP compatibilty


end
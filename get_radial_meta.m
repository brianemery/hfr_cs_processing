function R = get_radial_meta(R,APM,ftime,rkm)
% GET RADIAL META - meta data and hfrp struct formatting 
%
% Gets meta data and formats radial struct per HFRP formatting.
%
% Called by run_cs_processing*

% Copyright (C) 2019 Brian Emery

% CONVERT DOA STRUCT TO HFRP RADIAL STRUCT

R.Type = APM.Type(1:4);

R.SiteName = APM.SiteName;

R.SiteOrigin = APM.SiteOrigin;

R.SiteCode = 1;

% R.FileName = fileparts_name(fname);

R.TimeStamp = ftime; % fnames_to_times(fname,['CSS_' R.SiteName '_'],'yy_mm_dd_HHMM');

R.FileName = ['RDLm_' R.SiteName '_' datestr(R.TimeStamp,'yyyy_mm_dd_HHMM') '.mat'];

R.RangeBearHead = [];

% Generate Range 
R.RangeBearHead(:,1) = R.RngIdx(:,1) * rkm;


% CUSTOM FOR GLRT TESTING
R.Bear.MU = cwN2ccwE(R.Bear.MU); % <--- NOTE THIS IS NOT RANGEBEARHEAD!
R.Bear.ML = cwN2ccwE(R.Bear.ML);


% Idempodent fix I hope
R.LonLatUnits = R.LonLatUnits(1,:);
R.RangeBearHeadUnits = R.RangeBearHeadUnits(1,:);

% final cleanup - stuff not used on real data
R = rmfield(R,{'RmsBrgErr','BrgDiff','RomsBrg'});

% need this
R.ProcessingSteps{1,end+1} = 'get_radial_meta';



% DISABLE THE BELOW FOR GLRT TESTING ... which applies detection ...
%
% Code below is now being done by other things like apply_detection.m and
% by doa_to_radial_struct.m
return

% convert to bearing to ccwE
% Note from APM.README.BEAR_Units.BEAR_Units =  'degCWN'
R.RangeBearHead(:,2) = cwN2ccwE(R.Bear);

R.RangeBearHead(:,3) = NaN; % need this still 

% Generate LonLat
R.LonLat = rangeBear2LonLat(R.RangeBearHead,R.SiteOrigin);


R.RadComp = R.RadVel;



% % Populate U, V
% [R.U, R.V, R.Error, R.Flag] = deal(NaN(size(R.RadComp)));
% 
% % Need to compute heading ...
% % * HERE *
% 
% % From loadRDLfile, unchecked really 
% R.U = R.RadComp .* cosd(R.RangeBearHead(:,3));
% R.V = R.RadComp .* sind(R.RangeBearHead(:,3));

% % Add CS file name
% R.CSFileName = CS.FileName;


% final cleanup - dBear not used on real data
R = rmfield(R,{'Bear','RadVel','RngIdx','dBear'});


end


function R = get_radial_meta(R,APM,ftime,rkm)
% GET RADIAL META - meta data and hfrp struct formatting 
%
% Gets meta data and formats radial struct per HFRP formatting.
%
% Called by run_cs_processing*

% Copyright (C) 2019 Brian Emery

% 20220407 Adding code to handle empty struct input for R, which happens 

% CONVERT DOA STRUCT TO HFRP RADIAL STRUCT

R(1).Type = APM.Type(1:4);

R.SiteName = APM.SiteName;

R.SiteOrigin = APM.SiteOrigin;

R.SiteCode = 1;

% R.FileName = fileparts_name(fname);

R.TimeStamp = ftime; % fnames_to_times(fname,['CSS_' R.SiteName '_'],'yy_mm_dd_HHMM');

R.FileName = ['RDLm_' R.SiteName '_' datestr(R.TimeStamp,'yyyy_mm_dd_HHMM') '.mat'];

R.RangeBearHead = [];


if ~isempty(R.Bear)
    
    % Generate Range
    R.RangeBearHead(:,1) = R.RngIdx(:,1) * rkm;
    
    
    % CUSTOM FOR GLRT TESTING
    fn = fieldnames(R.Bear);
    for i = 1:numel(fn)
        R.Bear.(fn{i}) = cwN2ccwE(R.Bear.(fn{i})); % <--- NOTE THIS IS NOT RANGEBEARHEAD!
        % if isfield(R.Bear,'ML')
        %     R.Bear.ML = cwN2ccwE(R.Bear.ML);
    end
    
    
    % Idempodent fix I hope
    R.LonLatUnits = R.LonLatUnits(1,:);
    R.RangeBearHeadUnits = R.RangeBearHeadUnits(1,:);
    
end

% final cleanup - stuff not used on real data
R = rmfield(R,{'RmsBrgErr','BrgDiff','RomsBrg'});

% need this
R.ProcessingSteps{1,end+1} = 'get_radial_meta';



end


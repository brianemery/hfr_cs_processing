function R = get_radial_meta(R,APM,ftime,rkm)
% GET RADIAL META - meta data and hfrp struct formatting 
%
% Gets meta data and formats radial struct per HFRP formatting.
%
% Called by run_cs_processing*, but this is also used by Radar Simulation
% code on not-real data (eg check_plot_dales_stuff.m).
%
% Typically these run together:
% - get_radial_meta.m
% - apply_detection.m
% - doa_to_radial_struct.m or maybe detection_codar_post_proc.m


% Copyright (C) 2019 Brian Emery

% 20220407 Adding code to handle empty struct input for R, which happens 

% CONVERT DOA STRUCT TO HFRP RADIAL STRUCT

% special case
if nargin < 4, rkm =[]; end

R(1).Type = APM.Type(1:end);

R.SiteName = APM.SiteName;

R.SiteOrigin = APM.SiteOrigin;

R.SiteCode = 1;

% R.FileName = fileparts_name(fname);

R.TimeStamp = ftime; % fnames_to_times(fname,['CSS_' R.SiteName '_'],'yy_mm_dd_HHMM');

R.FileName = ['RDLm_' R.SiteName '_' datestr(R.TimeStamp,'yyyy_mm_dd_HHMM') '.mat'];

R.RangeBearHead = [];


if ~isempty(R.Bear)
    
    if ~isempty(rkm)
        % Generate Range
        R.RangeBearHead(:,1) = R.RngIdx(:,1) * rkm;
    else
        try R.RangeBearHead(:,1) = R.Range; catch, end
    end
    
    % CUSTOM FOR GLRT TESTING
    fn = fieldnames(R.Bear);
    for i = 1:numel(fn)
        R.Bear.(fn{i}) = cwN2ccwE(R.Bear.(fn{i})); % <--- NOTE THIS IS NOT RANGEBEARHEAD!
        % if isfield(R.Bear,'ML')
        %     R.Bear.ML = cwN2ccwE(R.Bear.ML);
    end
    
    
    % Idempodent fix I hope
    R.LonLatUnits = {'Decimal Degrees','Decimal Degrees'};
    R.RangeBearHeadUnits = {'km','Degrees_ccw_from_east','Degrees_ccw_from_east'};
    
end



% % final cleanup - stuff not used on real data
% % .. FUTURE maybe remove these if empty?
% R = rmfield(R,{'RmsBrgErr','BrgDiff','RomsBrg'});

% need this
R.ProcessingSteps{1,end+1} = 'get_radial_meta';



end


function radial_merge(wd,od)
% RADIAL MERGE - average subhourly radials to hourly like RadialMerger
% radial_merge(wd,od)
% 
% INPUT
% wd  - directory containing *.mat radials made from CSS files with
%       run_cs_processing.m 
% od  - output directory
%
% EXAMPLE
% 
% SEE ALSO
% bin_radials.m, bin_data_temporally.m, bin_data_struct.m

% Copyright (C) 2017 Brian Emery
%
% Version 12-Apr-2017 12:08:42


% STATUS: punt, see look_at_baseline_sni1_sci1.m


% check for test case
if strcmp('--t',wd), test_case, return, end


% get file list
flist = get_file_list(wd,'RDLm*');

% get site name
NM = cosFileNameParts(flist{1});

% get serial time array associated with flist
stime = fnames_to_times(flist,['RDLm_' NM.SiteName '_'],'yyyy_mm_dd_HHMM');



keyboard



end

function test_case
% TEST CASE
% 
% test case directory: /m_files/test_data/

wd = '/projects/error_covariance/data/sni1/RDLs/';
od = '/projects/error_covariance/data/sni1/RDLm/';

radial_merge(wd,od)


end

function t = fnames_to_times(flist,nmstr,dstr)
% FILE NAMES TO TIMES - more generic version of fname2time
% stime = fnames_to_times(flist,nmstr,dstr)
%
% INPUT
% flist  - cell list of file names, eg output of get_file_list.m
% nmstr  - part of file name to scrub off
% dstr   - date string format used in the name
%
% OUTPUT
% stime  - matlab serial time
%
% EXAMPLE
%
% flist =  { '/Users/codar/Desktop/tot_2015_05_19.mat'
%            '/Users/codar/Desktop/tot_2015_05_25.mat'};
% 
% stime = fnames_to_times(flist,'tot_','yyyy_mm_dd')
%
% ANOTHER EXAMPLE
%
% Custom code to match RDL_i*AGL1 for example (ie dont worry about RA)
%  nmstr = [RDLstr '.*' site '_']




% TO DO
% make this a replacement for fname2time or a function that it calls
%   char input support
%
% Custom code to match RDL_i*AGL1 for example (ie dont worry about RA)
%  nmstr = [RDLstr '.*' site '_']

t = datenum(regexprep(fileparts_name(flist),nmstr,''),dstr);


end
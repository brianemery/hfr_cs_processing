function totNames = times2totnames(varargin)
% TIMES2TOTNAMES.M - convert a time range to cell list of total file names
% totNames = times2totnames(st,en,drive);
% 
% Takes an input time array, or 2 end points, and converts
% them into names of total vector files which could then be loaded by a 
% different program. 'drive' is optional and can be used to specify 
% the directory path. The defaultis 'c:\data\totals\' on a pc, /data/totals
% otherwise.
%
% Used by mfiles such as get_total_data.m

% Copyright (C) 2013 Brian Emery 
% 
% originally 11mar03 


% note that the total velocity data base begins sept 1997!

if strcmp('--t',varargin{1}), test_case, return, end

% sort out the inputs
drive =[];
st = {};

for i = 1:size(varargin,2)

    if isnumeric(varargin{i}), st{i} = varargin{i}; end

    if ischar(varargin{i}), drive = varargin{i}; end

end

% convert st to numeric
st = [st{:}]';


% define data path if it's not given
if isempty(drive) & ispc
    
    drive='c:\data\totals\';
    
elseif isempty(drive) & ismac
    
    drive='/data/totals/';
end



% check for times before data availability
st(1) = max([st(1) datenum(1997,09,01)]);

% create the time array of days
st =floor(st(1)):floor(st(end));



% Handle crossing of time 2011_01_01


% OLD TOTAL FILES

% get the unique year-month strings
yrmo = unique(cellstr([datestr(st(st<datenum(2011,1,1)),11) datestr(st(st<datenum(2011,1,1)),5)]));

if ~isempty([yrmo{:}])
    % concat strings
    totNames = strcat({add_filesep(drive)},'tot_',yrmo,'.mat');
    
else
    totNames = [];
end
   


% NEW TOTAL FILES
yyyy_mm_dd = datestr(st(st>=datenum(2011,1,1)),'yyyy_mm_dd'); % keyboard


if ~isempty(yyyy_mm_dd)
    % concat strings and add to other totNames
    totNames = [totNames; strcat({add_filesep(drive)},'tot_',yyyy_mm_dd,'.mat')];
    
end

end

function test_case
% TEST CASE 
%
% in new format


% NEW TOTALS, SINGLE DATE INPUT

% define answer key:
key = {'C:\folder\on\a\pc/tot_2012_06_11.mat'};

% define and run the test
st = [datenum(2012,6,11,0,0,0):1/24:datenum(2012,6,11,23,0,0)];

drive = 'C:\folder\on\a\pc';

totNames = times2totnames(st,drive);

check_result(totNames, key, 'testing new totals')    



% OLD TOTALS, END POINT DATE INPUT

% define answer key:
key = {'/data/totals/tot_0011.mat'; '/data/totals/tot_0012.mat'; };

% define and run the test
st = [datenum(2000,11,1,0,0,0) datenum(2000,12,11,16,0,0)];

drive = '/data/totals/';

totNames = times2totnames(st,drive);

check_result(totNames, key, 'testing old totals')    



% CROSS TIME BETWEEN OLD AND NEW, AS 2 INPUTS

% define answer key:
key = {'C:\folder\on\a\pc/tot_1012.mat'; 'C:\folder\on\a\pc/tot_2011_01_01.mat'};


% Cross the time between old format and HFRP
st = datenum(2010,12,15);
en = datenum(2011,1,1,6,7,0);

drive = 'C:\folder\on\a\pc';

totNames = times2totnames(st,en,drive);

check_result(totNames, key, 'testing cross from old to new')    

end

function check_result(rslt, key, str)
% CHECK RESULTS 
%
% str is a description of the test

if isequal(rslt,key)
    disp([str ' ... ok'])
else
    disp([str ' ... NOT ok'])  , keyboard 
end


end




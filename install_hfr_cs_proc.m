% HFR CS PROCESSING - installation script
%
% INSTRUCTIONS
% - download and unzip the toolbox (eg via the clone or download button)
% - cd into the directory, eg:
%      cd ~/Downloads/hfr_cs_processing-master/
% - run this script:
%     >> install_hfr_cs_proc
%

% June 2020 Brian Emery

% get the directory locating this file
wd = fileparts( mfilename('fullpath') );

% list the directories
flist = dir(fullfile(wd, ['**' filesep '*.*']));

% get just the directories
flist = flist([flist.isdir]);  

% ... in a cell array
flist = unique({flist.folder}');

% scrub out .git and private
flist = flist(cellfun('isempty',regexp(flist,'.git')));
flist = flist(cellfun('isempty',regexp(flist,'private')));

% add these to the path
disp('Modifying MATLAB Path ...')

addpath(wd)

for i = 1:numel(flist)
    addpath(flist{i})
end

clear
function flist = get_file_list(wd,fstr)
% GET FILE LIST - get cell list of files including path
% flist = get_file_list(wd,fstr)
% 
% return empty cell if no files found
%
% wd can be a cell list of directories, but there is no need to 
% list the sub directories since this function handles the recursion
%
% EXAMPLE:
% flist = get_file_list(wd,'TRAK*');

% Copyright (C) 2013 Brian Emery

% Version 20130317 - added use of bash's find to do recursion

% Example of bash's find (with word count)
% find /Codar/SeaSonde/Archives/Spectra -name CSQ* | wc -l
%
% TO DO 
% Need to capture and display errors when I get the wd or fstr wrong
% Need to handle case where directory not found, get result = 'long
% string'
%
% Need to fix case where the find finds this, both of which fit: 
% flist(1:2)
% 
% ans = 
% 
%     '/projects/meas_vs_ideal/data/sni1/RDLm'
%     '/projects/meas_vs_ideal/data/sni1/RDLm/2014_03/RDLm_sni1_2014_03_13_0400.ruv'
%
% maybe test that it's a directory?
%
% this is a fix
% 
% % bug in getfile list
% tf = ones(numel(meas),1);
% 
% for i = 1:numel(meas)
%     tf(i) = exist(meas{i},'file');
% end
% 
% % make tf logical
% tf(tf == 2) = 1; tf(tf ~= 1) = 0; tf = logical(tf);
% 
% meas = meas(tf);




% check for test case
if strcmp('--t',wd), test_case, return, end


try
    % FAST WAY FOR MAC'S
    
    % check convert input to cell
    if ischar(wd)
        wd = cellstr(wd);
    end
    
    % start with empty flist
    flist = {};
    
    % loop over input directories
    for i = 1:numel(wd)
        
        % trival check for directory presence
        if isdir(wd{i})
            
            % the bash 'find' breaks when on the pwd, so change if needed
            if strcmp(pwd,wd{i}), [~,~] = system('cd ..'); end
            
            
            % remove trailing slash for clean output of 'find' below
            wd{i} = rm_slash(wd{i});
            
            % use bash's find to do recursion and generate list
            % sort is needed for linux file systems
            [status,result] = system(['find ' wd{i} ' -name ' '''' fstr ''' | sort']);
            
            
            % concat with previous, split results of find using eol
            if ~isempty(result) && status == 0
  
                flist = [flist; regexp(result(1:end-1),'\n','split')'];
                               
            end
            
            
        else
            
            disp([wd{i} ': Directory Not Found'])
            
        end
        
    end
            
catch E
    
    
    % OLD WAY - SHOULD WORK ON PC'S
    
    % check inputs, run recursive call if cell directory input
    if iscell(wd)
        flist =[];
        
        for i = 1:numel(wd)
            flist = [flist; get_file_list(wd{i},fstr)];
        end
        
        return
    end
    
    
    % add file separation
    wd = add_filesep(wd);
    
    % get list of files in order returned by O.S.
    flist = dir(strcat(wd,fstr));
    
    % output full path if files found
    if ~isempty(flist)
        flist = strcat(wd,cellstr(char(flist.name)));
    else
        flist = {};
    end
    
end



% check that all we've found are files

if ~isempty(flist)
    
    tf = ones(numel(flist),1);
    
    for i = 1:numel(flist)
        tf(i) = exist(flist{i},'file');
    end
    
    % make tf logical
    tf(tf == 2) = 1;
    tf(tf ~= 1) = 0;
    
    tf = logical(tf);
    
    % subset the list
    flist = flist(tf);
    
end





end

function wd = rm_slash(wd)
% wd is char

if strcmp(wd(end),filesep), wd = wd(1:end-1); end

end



function test_case
% TEST CASE

% flist base value (CELL)
base = {'/m_files/test_data/move_file/src/test_file.txt'};

% check listing works
flist = get_file_list('/m_files/test_data/move_file/src','test*');

check_result(flist, base, 'Test 1')    
    

% check empties returned if no files found
flist = get_file_list('/m_files/test_data/move_file/dest/','test*');

check_result(flist, {}, 'Test 2')


% check multiple input directories
wd ={'/m_files/test_data/move_file/src', '/m_files/test_data/move_file/dest',};
flist = get_file_list(wd,'test*');

check_result(flist, base, 'Test 3')


% test a good sized list for time ...
flist = get_file_list({'/m_files/test_data/cs_filter/good_cop_case', ...
                       '/m_files/test_data/cs_filter/bad_cop_case'},'CSQ*');

                   
% load previous result and check
load('/m_files/test_data/get_file_list.mat','prev')

check_result(flist, prev, 'Test 4')


% check empty returned for no directory found
flist = get_file_list('/m_files/test_data/move_file/dest/none','test*');

check_result(flist, {}, 'Test 5')


% need a check for recursion and strip off of dir names in the file list

keyboard
end

function check_result(flist, base, testStr)

if isequal(flist,base)
    disp([testStr ' passed'])
else
    disp([testStr ' failed'])  , keyboard 
end


end



% EXTRA CODE

function [fnames,fullfnames] = list_files(wdir,str,tf)
% LIST FILES - create cell list of directory contents
% [fnames,fullfnames] = list_files(wdir,str,tf)
%
% take cell input, optionally include directory ...
% with path included ...
%
% bug case:
%  [fnames,fullfnames] = list_files(wd,'STAT_sci1_2011_10_30_M_r*_d0_3x3.sdt')
% need to eliminate 2nd dir when no files are there ...
% 

% Set default input arguements
if nargin < 3, tf = true; end
if nargin < 2, str = '*'; end

% convert to cell inputs
if ischar(wdir)
    wdir = cellstr(wdir);
end

% init output cells
[fullfnames,fnames] = deal(cell(numel(wdir),1));

% GET FILE LIST
for i = 1:numel(wdir)
    
    nm = dir([wdir{i} str]); 
    
    % file names
    fnames{i} = {nm.name};
    
    % optionally append directory to file name
    if tf
        fullfnames{i} = strcat(wdir{i},fnames{i}); % strvcat(csqNames{:})
    end
end


% create outputs
fnames = [fnames{:}];
fullfnames = [fullfnames{:}];



end

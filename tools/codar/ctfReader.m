function S = ctfReader(file,tf)
% CTFREADER.M - Read Codar Table Format files (CTF) into a structure
%  S = ctfReader(file,swch);
%
% Uses column names given in header field 'TableColumnTypes' as the
% field names in the structure. Works on multiple files at a time
% (putting them into an indexed structure). Mimics some of the fields
% in a typical HFR_Progs stucture.
%
% Set swch to false (0) to prevent the automated saving of column format
% strings, eg for compiled code.
%
% See also ctfWriter.m, load_pattern_file.m, loadRDLFile.m, ctf_concat.m


% Copyright (C) 2011 Brian Emery
% Version 1.0, 25jun07
%         2.0  30jun08 with some borrowing from getRDLHeader.m
%         2.1  11Oct10 Clean up and details for compatibility with new
%                      ctfWriter.m. Now outputs footer.
%         2.2  18Jan12 added test case, detect precisions code
%              (ctf_formats.m and ctf_descriptions.m)

% TO DO
% - use [desc,units] = ctf_descriptions to get a LABS struct of plotting
%   labels
% - read addtional tables (using TableStart)

% NOTE
% TimeStamp is both a codar key name associated with the file creation
% time, and the name for the time array field in HFRProgs. Also, there is
% a footer name called ProcessedTimeStamp. From a .loop file for example:
%
% %TimeStamp: 2010 10 11  13 03 14
% %ProcessedTimeStamp: 2010 10 11  13 25 40
%
% meanwhile, some of the CTF files file names are a sort of time stamp.
% Try to sort this all out below


if strcmp('--t',file), test_case, return, end


S=struct();

if ischar(file)
    file = cellstr(file);
end

% % check file is found
% if ~(exist(file,'file') == 2)
%     disp(['ctfReader.m: file not found (' file ')'])
%     S = []; return
% end
%



% loop over file names
for i = 1:numel(file)
    
    % get some info from the file name
    NM = cosFileNameParts(file{i});
    
    [pth,fname,ext] = fileparts(file{i});
    
    % Add standard fields to top of structure
    S(i).Type = NM.Type;
    S(i).SiteName = NM.SiteName;
    S(i).SiteOrigin = [];
    S(i).CreateTimeStamp = datestr(NM.TimeStamp);
    S(i).StructTimeStamp = creation_info;
    S(i).FileName = [fname ext];
    
    
    
    % READ DATA
    % experimental code here for dealing with bad line encodings ...
    
    try
        % [^...] - reads characters not matching characters between the
        % brackets until first matching character or whitespace
        % (returns cellstr)
        dat = textread(file{i},'%[^\n]',-1);
        
        
    catch E
        
        disp(E.identifier)
        disp('retry ...')
        
        % try after applying unix character translate tool
        [stat,rslt] = system(['tr -d ''\000'' < ' file{i} ' > ' file{i} '.tr']);
        
        dat = textread([file{i} '.tr'],'%[^\n]',-1);
        
    end
    
    % Get the header and footer
    % find text fields
    isTxt = strncmp('%',dat,1);
    idx   = [1:numel(dat)]';
    
    % Find footer if it exists
    ft = strmatch('%TableEnd:',dat);
    if isempty(ft), ft = numel(dat); end
    
    % Apply indexing
    S(i).Header = dat(isTxt & idx < ft(1));
    S(i).Footer = dat(isTxt & idx >= ft(end));
    
    
    % Kaplan's code to get variable key words 'names' and 'values' from header
    % Might be better to do S(i).(HdrName) = HdrValue?
    [S(i).HdrNames,S(i).HdrValues] = get_names_values(S(i).Header);
    
    
    
    % PUT DATA IN STRUCTURE
    
    % create empty time array now to keep things tidy
    S(i).TimeStamp = [];
    
    % Put remaining data into a matrix
    % this is slow but robust
    % dataStr = char(dat(~isTxt));
    % dat = str2num( dataStr );
    %
    % ... and this is about 6x faster and also robust (gets same thing),
    % but we need to keep a few of the lines for ctf_formats
    ix = find(~isTxt);
    
    % load the data if there is any
    if ~isempty(ix)
        
        try
            dataStr = char(dat(ix(1:10)));
            dat = load(file{i});
            
        catch E
            % catch badly written files, specifically, ones that are quit while
            % writting the last line
            % ... maybe this? A = textscan(fid, ,'Headerlines', min(ix)-1)
            disp(E.message)
            disp('... trying other method ...')
            % dat = str2num( char( dat(ix(1:end-1)) ) );
            
            % get number of columns ... all in try catch?
            [~,~,colstr] = getNameValuePair( 'TableColumns', S.HdrNames, S.HdrValues);
            
            % check ismac first?
            
            % get unique temp file name
            %tmp = [prefdir '/tmp_' datestr(now,30) '.txt'];
            tmp = ['/private/tmp/tmp_' datestr(now,30) '.txt'];
            
            % Create Tmp file
            % run awk command, eg: awk 'NF>=29' LOOP_ptc1_19001122_000000.loop
            % direct output to tmp file in pwd?
            [s,r] = system(['awk ''NF>=' colstr{1} ''' ' file{1} ' > ' tmp ]);
            
            % now try load
            dat = load(tmp);
            
            % remove tmp file
            [r,s] = system(['rm ' tmp]);
            
        end
        
        % Use field names from 'TableColumnTypes'. Note that it may occur more than
        % once, which is the case with radial files.
        % example of val = 'BEAR A13R A13I A23R A23I FLAG A13M A13P A23M A23P'
        [ix,nm,val]=getNameValuePair('TableColumnTypes',S(i).HdrNames,S(i).HdrValues);
        fn = cellstr(strparser(val{1}));
        
        for ii = 1:numel(fn)
            
            % Catch error, put in empty field and if data not found
            try
                S(i).(fn{ii}) = dat(:,ii);
            catch
                S(i).(fn{ii}) = [];
            end
            
        end
        
        
        % ADDITIONAL META DATA
        
        % Add meta info to bottom of struct
        S(i).ProcessingSteps = {mfilename};
        
        % Try to create a time array for timeseries CTF data
        try
            S(i).TimeStamp = datenum(S(i).TYRS,S(i).TMON,S(i).TDAY,S(i).THRS,S(i).TMIN,S(i).TSEC);
        catch E, end
        
        % Try to use the 'TimeStamp' from the CTF file as the creation date
        try
            S(i).CreateTimeStamp = S(i).HdrValues{strcmp('TimeStamp',S(i).HdrNames)};
        catch E, end
        
        % Try to get the site location
        try
            S(i).SiteOrigin([2 1]) = str2num(S(i).HdrValues{strcmp('ReceiverLocation',S(i).HdrNames)});
        catch E, end
        
        
        % Try to store precisions for ctfWriter, but dont allow errors to stop
        % the execution
        if nargin < 2 || (nargin > 1 && tf)
            try
                ctf_formats( dataStr, fn);
                disp(['Saved CTF 4 char key data to database'])
                
            catch E
                % warning('ctfReader.m call to ctf_formats.m failure')
            end
        end
        
        
        % Try to store column descriptions for ctfWriter also
        % this is not going to work at this time with STAT files!
        if nargin < 2 || (nargin > 1 && tf)
            try
                [S(i).LABS,S(i).Units] = ctf_descriptions( S(i).HdrNames, S(i).HdrValues);
            catch E
                % warning('ctfReader.m call to ctf_descriptions.m failure')
                % keyboard
            end
        end
    end
end
end

function [names,values] = get_names_values(hdr)
% GET NAMES VALUES - get name and value pairs from CTF format header
%
% EXAMPLE
% Some name, value pairs:
% %ReceiverLocation:  34.4078333 -119.8783333
% %RangeCells: 88
% %DopplerCells: 32
% %DopplerUpdateRate: 4
%
%
% From David Kaplan's code to get variable 'names' and 'values' from header
% Next split lines if desired.
%
% from getNameValuePair.m

[names, values] = deal(cell(length(hdr),1));

for k = 1:length(hdr)
    % This way is considerably more efficient than using strtok
    ii = min( [ find( hdr{k} == ':' ), length(hdr{k})+1 ] );
    names{k,1} = strtrim(hdr{k}(2:ii-1)); % Remove initial %
    values{k,1} = strtrim(hdr{k}(ii+1:end)); % Removie initial :
    
end

end

function test_case
% TEST CASE
%

% TEST BAD CHARACTER
wd = '/m_files/test_data/ctfReader/';
fn = 'LOOP_ARG1_130418_000000.loop';

T = ctfReader([wd fn]);

% clean up
[r,s] = system(['rm ' wd fn '.tr']);

keyboard

% TEST IRREGULAR COLUMNS
wd = '/m_files/test_data/ctfReader/';
fn = 'LOOP_ptc1_19001122_000000.loop';

T = ctfReader([wd fn]);

% should be the same as this file:
S = ctfReader([wd 'LOOP_ptc1_19001122_000001.loop';]);

test_check(S,T)

keyboard

% BASIC TEST

T = ctfReader(['/m_files/test_data/ctfWriter/LOOP_cop1_100818_192848.loop']);

% compare with previous (verified) output
load /m_files/test_data/ctfReader/LOOP_cop1_100818_192848.mat


test_check(S,T)


keyboard


% USE THIS LOOP FILE TO POPULATE CTF_FORMATS
T = ctfReader(['/m_files/test_data/ctfReader/LOOP_new1_000000_000000.loop']);
T = ctfReader(['/m_files/test_data/ctfReader/LOOP_sni1_20120223_000000.loop']);


end

function test_check(S,T)

% remove time and other independent content
fn = {'StructTimeStamp','FileName'};
T = rmfield(T,fn);
S = rmfield(S,fn);


if isequalwithequalnans(S,T)
    disp(['ctfReader.m: test case passed'])
    
else
    disp(['ctfReader.m: test case NOT passed'])
    
    keyboard
    
    fn = fieldnames(S);
    
    for i = 1:numel(fn)
        
        if ~isequal(S.(fn{i}),T.(fn{i}))
            keyboard,
        end
        
    end
    
    keyboard
    
end

end


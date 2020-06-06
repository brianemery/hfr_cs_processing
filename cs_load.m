function CSL = cs_load(flist,ftimes,rtime,tau,CSL)
% CS LOAD - load CS data for radial processing
% 
% CSL = cs_load(flist,ftimes,rtime,tau,CS,CSL)
%
% Called by run_cs_processing.m, this manages the rolling collection of CS
% data to average togther, for example to make CSS from CSQ. It outputs
% the multi-element struct of individual CS observations.
% 
% INPUTS
% flist  - full list of CS files to load and or average
% ftimes - times of files in flist
% rtime  - time of radial file output associated with the +/- tau CS files
% tau    - averaging period in minuts
% CSL    - struct of everything in the average
% 
% Set tau = 0 for no averaging - just load one (TO DO)
%
%
% REFERENCE
% http://www.analog.com/media/en/technical-documentation/dsp-book/dsp_book_Ch15.pdf

% Copyright(C) 2017 Brian Emery


% need to make pass through for CSS ?
if tau == 0, keyboard, end


% get list of files already in the struct
inlist = {CSL.FileName}';

% find times that should be in the average
cx = find( ftimes >= ( rtime - tau/(2*1440) ) & ftimes <= ( rtime + tau/(2*1440) ) );

% get list of files to remove from the CSL struct
rmlist = setdiff(inlist,flist(cx));


% find files to add to the list
addlist = setdiff(flist(cx),inlist);


% load the new data 
if ~isempty(addlist)

    for i = 1:numel(addlist)
        % CSL(end+1) = ReadCS(addlist{i});
        CSL(end+1) = load_by_type(addlist{i});
    end
    
end


% remove data not needed for the average over tau
CSL = CSL( ~ismember([inlist; addlist],rmlist) );



return

% OLDER VERSION

% Ok ... why not just keep the rolling struct (CSL) and recompute the
% average each time? should be fast enough ... 




% get the substruct to remove
CSrm = CSL( ismember(inlist,rmlist) );
    
% load the new data to add
if isempty(addlist)
    CSn = struct([]);
else
    CSn = cs_struct(numel(addlist));
    for i = 1:numel(addlist)
        CSn(i) = ReadCS(addlist{i});
    end
end

% compute the new moving average
CS = cs_mv_average(CS,CSn,CSrm);

% get rid of the old 
CSL = CSL( ~ismember(inlist,rmlist) );

% add the new
if ~isempty(CSn)
    CSL(end+1:end+numel(CSn)) = CSn;
end

% % Check plotting
% figure(1), hold on, plot(1:512,10*log10(CS.antenna3Self(:,10)),'-')
%     
% keyboard    
% 
% % at any time the average of CSL should be the same as CS
% CSA = cs_average(CSL);
% plot(1:512,10*log10(CSA.antenna3Self(:,10)),'-b')



end

function CS = load_by_type(fname)

[~,~,ext] = fileparts(fname);

if strcmp('.cs',ext) 
    
    % S = ReadCS(fname);
    CS = cs_read(fname);
    
elseif strcmp('.mat',ext) 
    
    CS = load(fname,'CS'); CS = CS.CS; %jfc
end




end

function CS = cs_mv_average(CS,CSn,CSrm)
% CS MV AVERAGE 
%
% local code for maintaining a running average. 
%
% Breaks if CSn is empty, need to detect and handle appropriately
% also if CSrm empty ... 1578
% ... this handled outside which sets structs to empty, they then get numel
% = 0, and that fixes the addition stuff


% get number in the average
n = CS.N;

% field names to act on
fn = {'antenna1Self', 'antenna2Self', 'antenna3Self', ...
       'antenna12CrossSp','antenna13CrossSp','antenna23CrossSp'};
   
% number to remove 
m = numel(CSrm);

% number to add
p = numel(CSn);


for i = 1:numel(fn)
    
    % multiply
    CS.(fn{i}) = CS.(fn{i}).*n;
    
    % subtract
    for j = 1:m
        CS.(fn{i}) = CS.(fn{i}) - CSrm(j).(fn{i});
    end
    
    % add
    for k = 1:p
        CS.(fn{i}) = CS.(fn{i}) + CSn(k).(fn{i});
    end
    
    % divide by the new N
    CS.(fn{i}) = CS.(fn{i})./(n-m+p);
end
   
% Update output
CS.N = n-m+p;

% Use newest header also
try
    CS.Header = CSn(end).Header;
catch 
end


end 


% OBSOLETE I THINK
function get_averaging_setup
% % CSQ Averaging ...
% % Should break this out for testing and stuff
% function [ftimes,rtimes,CS,CSL] = get_averaging_setup(site,flist,tau)
% GET AVERAGING SETUP - get file times and load inital files 
%
% Converts file list of names to file times
% Loads initial 
%
% INPUTS
% site  -  char site name
% flist -  list of cs files
% tau   -  averaging time in minutes
%
% OUTPUTS
% ftimes - file times from the file names
% rtimes - times to output radial files
% CS     - average of N CS data files over tau
% CSL    - N-element struct the N CS data files
  
keyboard

% assume using csq's if averaging
 
if tau < 1
    
    % CSS case --- need to test/develop 
    ftimes = sort(fnames_to_times(flist,['CSS_' site '_'],'yy_mm_dd_HHMMSS'));
    rtimes = ftimes;
    [CS,CSL] = deal([]);
    
else
    
    % CSQ case
    ftimes = sort(fnames_to_times(flist,['CSQ_' site '_'],'yy_mm_dd_HHMMSS'));
         
    % find ones to load
    cx = find( ftimes <= ( ftimes(1) + tau/1440 ) );
        
    % preallocate
    n = length(cx);
    
    CSL = cs_struct(n);
    
    % load first batch
    for i = 1:n
        CSL(i) = ReadCS(flist{cx(i)});
    end
            
    % average
    CS = cs_average(CSL,n); CS.N = n;
           
    % get times to output radial files (every 10 min)
    rtimes = unique(round(ftimes(n/2:end)*144)/144);
    
end

end


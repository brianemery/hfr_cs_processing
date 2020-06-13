function run_cs_processing %(site)
% RUN CS PROCESSING - cross spectra to radials driver script
%
% Typically, this is followed by run_cs_post_processing.m, which applies
% the detection algorithms, updates meta data and applies some formatting.
%
% ARCHITECTURE IDEAS
%
% - could follow range processing
% - arbitrary array geometry, fft length, etc
% - arbitrary doa method 
% - future use of things like different detectors, MAP algorithm, ...
% - tests using Anthony's data for validation, but also could build some
%   edge case tests using simiulations (eg 0-360 transition, etc)
% - simplify for users: standard directory structure? Separate config
%   creation (CFG struct as input).
% 
%
% DEFAULT BEHAVIOR
% - Outputs radials every 10 min
% - Uses the list of CS files to set processing times

% Copyright (C) 2017-2020 Brian Emery
%
% Last update 28-Feb-2020 11:40:23


% TO DO
% - Store meta data in the data structure ! all this config stuff !!!
% - guess I need to specify a header.txt file to get the range increment
%   since it's not generally available anywhere else. For now assume that
%   distance to first range cell is same as the spacing (given in rkm)
% - add some user feedback "processing ... blah"
% - Document processing settings: Tau, FOL method and settings, site,
%   flist, APM used, what else ...
% 
% - Need test data set ... one set of CSQ and R to compare new with prev
%   maybe compare with a codar version of the radial ...
%
% - I have a 100cm/s cutoff, but sometimes get larger Vr - why?
%
% CHECK FOR RNG INDEX ERRORS
% - some code assumes continuous numbers of range cells ... need to not
% assume that eg the FOL cell arrays are indexed same as the number of
% range cells


% NOTE: start parpool or matlabpool (see line 42 of doa_on_cs.m)... on
% knot?
%
% NOTE2: sci1 on vici with the below settings is the defacto test ... need
% to encode it ...
%
% NOTE3: Anthony uses 1024 FFTs with 50% ... or maybe 2048 FFTs with 78%
% overlap, in 30 min 


% CONFIGURE

site = 'cop1';  % same as site file name in radial file
tau = 27;       % Set the averaging time in minutes 

% Main Data Folder (Has CS, RDL subdirs and APMs in it)  
%wd = ['/scratch_local/emery/' site]; % on Knot
% wd = ['/home/emery/data/' site];
% wd = ['/projects/detection/data/' site ];  
wd = '/projects/hf_winds/data/cop1';

% apm file to use ... need a better way to organize these ...
%apm_file = [wd '/MeasPattern_sci1_20130815.txt']; %'/MeasPattern_sni1_20130831.txt']; % 
% apm_file = [wd '/MeasPattern_ptc1_20170603.txt']; % /MeasPattern.txt']; %
apm_file = [wd '/MeasPattern_cop1_20181030.txt'];

% Define CS file directory
csd = [wd '/csq/'];

% % K = 3 4 6 10
% tau = [13 18 27 43];
% if nargin < 1, tau = 27; end

% Define radial data output directory ... gets made by code below
% rwd = [wd '/RDLs_K-' num2str(round(tau*60/256)) '/'];
rwd = [wd '/RDLs_zero_doppler/'];

% max number of emitters to search for
Nemit = 2; % probably 5 for other arrays


% FOL user parameters (see imageFOL.m help)
user_param = [40 100 10]; % [stdev maxcurr snrcutoff] <-- maybe use 10 for SNR? 

% Determine time (min) of CS to average, set to zero if no averaging (eg if
% processing CSS). Note that I use a sub-function to estimate K based on the
% number of actual input files, assuming CSQs
cs_type = 'CSQ';




% PRE - PROCESS

% un pack the struct for backward compatibility
% struct_unpack(CFG)

% get the apm
APM = load_pattern_file(apm_file);

if strcmp(site,'sni1')
    APM.SiteName = site;
    
    % For SNI, fix the 360-0 transition ... need a general solution to this
    APM.BEAR(APM.BEAR>180) = APM.BEAR(APM.BEAR>180)-360;
end

% list cs files (csq or css)
flist = get_file_list(csd,'CS*');


% convert file list to their times
ftimes = sort(fnames_to_times(flist,[cs_type '_' site '_'],'yy_mm_dd_HHMMSS'));
 
% get times to output radial files (this is every 10 min)
rtimes = unique(round(ftimes*144)/144); % rtimes = ftimes for CSS

   
% Get setup data
% get imageFOL inputs and info from CS file
[iBragg,v_incr,rkm] = get_setup_data(flist{1}); %<-- CHECK VS HEADER.TxT

% Check for and or make output directory
[~,~,~] = mkdir(rwd);

% init empty CS structs 
CSL = cs_struct(0);



% diff with previous list

% check save directory is there
dir_check(rwd)


% PROCESSING LOOP


for i = 1:numel(rtimes)  
    
    % load cs files based on what is needed to average together
    CSL = cs_load(flist,ftimes,rtimes(i),tau,CSL);    
    
        
    % ship removal step ...
    % this needs to happen before computing the average
    % needs min amount of cs data ... so punt
    try
        CSL = cs_ship_rm(CSL);
    catch E
        disp('cs_ship_rm error ...')
    end
    
    
    
    % compute the average of the CS for processing. 'K' is actually the
    % number of files, but a rough est of the snapshots
    [CS,K] = cs_average(CSL); 
        
    % add standard CS info 
    [CS.freqs,CS.Vrad,CS.fb] = getVelocities(CS.Header);
    
    % do this now, use the data in the doa subfunction
    CS.SNR = get_SNR(CS);

    
    disp('RUNNING ZERO DOPPLER PROCESSING')
    % % FOL processing (rows are range cells)
    % [FOreg, FOregi, Alims] = imageFOLs(CS.antenna3Self.',iBragg,v_incr,user_param);
    %
    %
    % % convert FOregi to peakIdx
    % peakIdx = convert_fol_ix(FOregi);
    
    % Holly and Hondo ... 256 is zero Doppler velocity, Holly is just into
    % RC 3, and Hondo is just between RC 14,15
    peakIdx([2 3 14 15],:) = {256} ;
    
    %  pcolor(10*log10(CS.antenna1Self.'))
    %  shading flat; colorbar; caxis([-160 -80]);
    %  hold on
    % h= plot(FOregi(:,2),FOregi(:,1),'c.');
    %
    
    % Estimate K based on data that has been loaded
    % K = estimate_k(CSL); % Now done by cs_average.m
   
    
    % make sure peakIdx is a cell array (each range cell a cell}    
    % Run DOA processing
    % (music or mle, music metrics, music error, etc)
    S = doa_on_cs(CS,APM,peakIdx,K,Nemit);  
    
              
    % PREVIOUSLY ...
%     % detection step - currently uses codar's music parameters
%     % ML = apply_test_result(ML); % try to skip this for ML ...
%     MU = apply_test_result(MU);
%     
%     % explicitly setting this here
%     MU.ProcessingSteps = {'imageFOLs.m','music.m','run_param_test.m'};
%     ML.ProcessingSteps = {'imageFOLs.m','mle_ap.m'}; %,'run_param_test.m'};
%     
%     
%     % GET RADIAL STRUCT OUTPUT
%     % clean up radial struct meta data
%     % Rmu = get_radial_meta(MU,APM,rtimes(i),rkm);
%     %Rml = get_radial_meta(ML,APM,rtimes(i),rkm);
%     
%    
%     % write the radial struct to a mat file
%     save([rwd Rmu.FileName],'Rmu','ML','APM','rtimes','rkm','i','K') %'Rml')


    % SAVE PRE-APPLICATION OF DETECTION
    % ... currently set to change brg to ccwE
    S = get_radial_meta(S,APM,rtimes(i),rkm);
    
    % store this - same for whole CS file
    S.K = K;
    
    save([rwd S.FileName],'S','APM','rtimes','rkm','i','K') 

end


return

keyboard
% ... can to post process on the ML data without detection


% POST PROCESS

% Concat structs ... make the temporal averages
% gather_radial_data(wd,site) ... save different files for the different
% DOA methods
run_cs_post_processing(wd,site)



disp('... done')

keyboard



end


function [iBragg,v_incr,rkm] = get_setup_data(fname)
% GET FOL INPUTS
%
% get stuff needed for imageFOL: Index of Bragg lines, and radial ocean
% current velocity increment in m/s

[~,~,ext] = fileparts(fname);

if strcmp('.cs',ext) 
    %CS = ReadCS(fname);
     CS = cs_read(fname);
elseif strcmp('.mat',ext) 
    CS = load(fname,'CS'); CS = CS.CS;
else
    keyboard
end

[~,Vrad,~] = getVelocities(CS.Header);

[~,~,dv] = getDopplerVelocities(CS.Header);

n = length(Vrad)/2;


% get approx indecies of Bragg lines
[~,iBragg(1)] = min( abs( Vrad(1:n)) ) ;
[~,iBragg(2)] = min( abs( Vrad(n:end)) ) ; 

% adjust last result
iBragg(2) = iBragg(2) + n - 1;

% velocity increment m/s
v_incr = dv/100; %mode(diff(Vrad))/100;

% get distance to first range cell
rkm = CS.Header.distToFirstRangeCell;

end

function K = estimate_k(CSL) % THIS IS NOW OBSOLETE see cs_average.m
% ESTIMATE K - get number of snapshots
%
% assumes files in CSL are CSQ files ...
%
% Compute or specify K, the number of independent samples used to form the
% covariance matrix ... this is used by music_error.m. According to
% ISS_CH4_SpectraToRadials.pdf, CSQs do not overlap.
%
% Estimate K based on the number of actual input files, rather than based
% on the times ...

% If CSL has an empty field, the cat will route around it and this will get
% the true number of data files, which would be different from numel(CSL)

K = size( cat(3,CSL.antenna3Self) ,3);



end

function dir_check(wdir)
% qad, from saveAsMakeDirectory ...

if ~isdir(wdir)
    
    % make the directory
    [s,~,~] = mkdir(wdir);
    
    if s, disp(['Created ' wdir]), end
end



end


function previous_configuations
% PUT THIS HERE FOR REUSE LATER ...

%site = 'ptc1'; %'sci1'; %'sni1'; %  same as site file name in radial file
site = tau; % repurpose tau for use as site string
tau = 27;

% Main Data Folder (Has CS, RDL subdirs and APMs in it)  
%wd = ['/scratch_local/emery/' site]; % on Knot
% wd = ['/home/emery/data/' site];
wd = ['/projects/detection/data/' site ];  

% apm file to use ... need a better way to organize these ...
%apm_file = [wd '/MeasPattern_sci1_20130815.txt']; %'/MeasPattern_sni1_20130831.txt']; % 
apm_file = [wd '/MeasPattern.txt']; %_ptc1_20170603.txt'];


% Define CS file directory
csd = [wd '/csq/'];

% % K = 3 4 6 10
% tau = [13 18 27 43];
% if nargin < 1, tau = 27; end

% Define radial data output directory ... gets made by code below
rwd = [wd '/RDLs_K-' num2str(round(tau*60/256)) '/'];

% max number of emitters to search for
Nemit = 2; % probably 5 for other arrays


% FOL user parameters (see imageFOL.m help)
user_param = [40 100 5]; % [stdev maxcurr snrcutoff]

% Determine time (min) of CS to average, set to zero if no averaging (eg if
% processing CSS). Note that I use a sub-function to estimate K based on the
% number of actual input files, assuming CSQs
cs_type = 'CSQ';





end
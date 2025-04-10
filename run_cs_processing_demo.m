function run_cs_processing_demo
% RUN CS PROCESSING DEMO - cross spectra to radials driver script demo
%
% Demonstrates a typcical usage of the toolbox and also acts as a test to
% verify that it works 
%
% USAGE
% set the working directory (wd) below to the directory containing data
% from BML1 included with the toolbox
% 
% If you have the parallel processing toolbox you can edit doa_on_cs.m
% (near line 80) to use it.
%
% The configuration settings below are typical for SeaSondes. 

% Copyright (C) 2020 Brian Emery



% CONFIGURE

% get the directory locating this file to use as the root directory
% (Has CS, RDL subdirs and APMs in it)
wd = fileparts( mfilename('fullpath') );

site = 'BML1';  % same as site file name in radial file

wd = [wd '/data/' site ];


% apm file to use
apm_file = [wd '/APM/MeasPattern_BML1.txt'];

% Define CS file directory
csd = [wd '/CSS/'];


% Define radial data output directory ... gets made by code below
rwd = [wd '/RDLs'];

% max number of emitters to search for
Nemit = 2; % for Seaonde; Probably 5 max for other arrays


% FOL (first order line) user parameters (see imageFOL.m help)
user_param = [40 100 10]; % [stdev maxcurr snrcutoff]  use 10 for SNR?

% Determine time (min) of CS to average, set to zero if no averaging (eg if
% processing CSS). Note that I use a sub-function to estimate snapshots (K)
% based on the number of actual input files, when using CSQs
tau = 30;
cs_type = 'CSS';

% if using 256 second FFTs:
% tau = [13 18 27 43];
% K = 3 4 6 10
%
% For this demo, set K = 7 (number of independent FFT's approx)
% ... important only for the music error calculation, and the likelihood
% ratio calculation
K = 7;

% Detection settings
% SeaSonde MUSIC parameters are set to [40 20 2] in doa_on_range_cell.m
% Set LR threshold here (between 10 and 40 for SeaSondes)
cut = 30;




% PRE - PROCESS

% un pack the struct for backward compatibility
% struct_unpack(CFG)

% get the apm - requires first 4 characters in the name to be 'Meas',
% 'Idea', 'SEAS', or 'LOOP'.
APM = load_pattern_file(apm_file);


% list cs files (csq or css)
flist = get_file_list(csd,'CS*');


% convert file list to their times
% For CSQ use this:
% ftimes = sort(fnames_to_times(flist,[cs_type '_' site '_'],'yy_mm_dd_HHMMSS'));
% For CSS use this:
ftimes = sort(fnames_to_times(flist,[cs_type '_' site '_'],'yy_mm_dd_HHMM'));

% Define the times to generate output radial files
% (this is every 30 min):
rtimes = ftimes(1): 30/1440 : ftimes(end);

% get imageFOL inputs and info from CS file
[iBragg,v_incr,rkm] = get_setup_data(flist{1}); %<-- CHECK VS HEADER.TxT

% Check for and or make output directory
[~,~,~] = mkdir(rwd);

% init empty CS structs
CSL = cs_struct(0);

% Gather some meta data to save later
CFG.README = 'Meta Data for run_cs_processing_demo.m';
CFG = struct_pack({'apm_file', 'user_param','tau','cs_type','K','cut','Nemit'},CFG);

% set other defaults
CFG.dmth = {'mu','ml'};
CFG.use_parfor = false;
CFG.mus_param = [10 5 8];




% check save directory is there
dir_check(rwd)


% PROCESSING LOOP


for i = 1:numel(rtimes)
    
    % load cs files based on what is needed to average together
    CSL = cs_load(flist,ftimes,rtimes(i),tau,CSL);
    
    % if CSL is empty ... skip to next i
    if ~isempty(CSL)
        
        
        % ship removal step ...
        % this needs to happen before computing the average
        % needs min amount of cs data ... so punt
        try
            CSL = cs_ship_rm(CSL);
        catch E
            disp('cs_ship_rm warning ...')
        end
        
        
        
        % compute the average of the CS for processing.
        [CS,~] = cs_average(CSL);
        
        % add standard CS info
        [CS.freqs,CS.Vrad,CS.fb] = getVelocities(CS.Header);
        
        % do this now, use the data in the doa subfunction
        CS.SNR = get_SNR(CS);
        
        
        % FOL processing (rows are range cells)
        [FOreg, FOregi, Alims] = imageFOLs(CS.antenna3Self.',iBragg,v_incr,user_param);
        
        
        % convert FOregi to peakIdx
        peakIdx = convert_fol_ix(FOregi);
        
        
        % make sure peakIdx is a cell array (each range cell a cell}
        % Run DOA processing
        % (music or mle, music metrics, music error, etc)
        S = doa_on_cs(CS,APM,peakIdx,K,CFG);
        
        % fill in some meta data, including Range, file name, etc
        S = get_radial_meta(S,APM,rtimes(i),rkm);
        
        
        
        
        % DETECTION STEP
        
        % use codar's music parameters
        
        D = apply_detection(S,S.Dual+1,'seasonde');
        
        D = doa_to_radial_struct(D,'MU');
        
        D = fix_multicol_fields(D);
        
        U = D; clear D
        
        
        % likelihood ratio detection
        % apply detection, seperate paths for Music and MLE
        LR = -2*log(S.LR.ML);  % exp(30/-2) ans = 3.0590e-07
        
        % compute emitters
        em = compute_emitters_from_lr(LR,cut);
        
        % apply detection
        D = apply_detection(S,em,'glrt'); %,tf);
        
        % house keeping
        D = doa_to_radial_struct(D,'ML');
        
        D = fix_multicol_fields(D);
        
        % store and clear
        L = D; clear D
        
        
        % MORE CLEAN UP
        
        % explicitly setting these here
        U.ProcessingSteps = {'imageFOLs.m','music.m','run_param_test.m'};
        L.ProcessingSteps = {'imageFOLs.m','mle_ap.m','glrt.m'};
        U.K = K; L.K = K;
        
        U.README = rmfield(U.README,{'Bear','RunTime','LR','Rm','Idx'});
        U = rmfield(U,{'RngIdx','RunTime','LR','Rm','Idx'});
        
        L = rmfield(L,{'RngIdx','RunTime','LR','Rm','Idx'});
        L.README = rmfield(L.README,{'Bear','RunTime','LR','Rm','Idx'});
        
        
        
        % write the radial struct to a mat file
        save([rwd '/' U.FileName{1}],'L','U','CFG')
        
        clear L U S
        
    end
    
end

disp('... done')


% Compare with the SeaSonde Processing file

R = loadRDLFile([wd '/RDLm/RDLm_BML1_2019_02_17_2100.ruv']);

ix = find(R.RangeBearHead(:,1) > 10 & R.RangeBearHead(:,1) < 12);

plot(R.RangeBearHead(ix,2),R.RadComp(ix),'-o')

hold on

S = load([wd '/RDLs/RDLm_BML1_2019_02_17_2100.mat']);

% look at music processed ...
iy = find(S.U.RangeBearHead(:,1) > 10 & S.U.RangeBearHead(:,1) < 12);
plot(S.U.RangeBearHead(iy,2),S.U.RadComp(iy),'*')

% ... with bearing uncertainties
errorbarx(S.U.RangeBearHead(iy,2),S.U.RadComp(iy),S.U.Err(iy));

% MLE processed
iz = find(S.L.RangeBearHead(:,1) > 10 & S.L.RangeBearHead(:,1) < 12);
plot(S.L.RangeBearHead(iz,2),S.L.RadComp(iz),'^')

xlabel('Bearing')
ylabel('Radial Velocity')

legend('SeaSonde Std Proc','HFR CS Proc - MUSIC', ...
       'HFR CS Proc - MUSIC Uncert','HFR CS Proc - MLE')


% look at some spectra

 flist = get_file_list([wd '/CSS/'],'CSS*');

for i = 19 %:25
    CS = cs_read(flist{i}); s= cs_plot_map(CS); %pause   
end

keyboard


end


function [iBragg,v_incr,rkm] = get_setup_data(fname)
% GET FOL INPUTS
%
% get stuff needed for imageFOL: Index of Bragg lines, and radial ocean
% current velocity increment in m/s

[~,~,ext] = fileparts(fname);

if strcmp('.cs',ext)
    % CS = ReadCS(fname);
    CS = cs_read(fname); % CS Ver 6 compatible
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

function dir_check(wdir)
% qad, from saveAsMakeDirectory ...

if ~isdir(wdir)
    
    % make the directory
    [s,~,~] = mkdir(wdir);
    
    if s, disp(['Created ' wdir]), end
end



end

function  D = doa_to_radial_struct(D,doa_field)
% DOA TO RADIAL STRUCT - final clean up at this point, but this could be a
% function that does the whole thing start to end ...
%
% doa_field = 'ML' or 'MU'

% keep these around for later?
D.OtherMatrixVars.eigValues = D.eigValues;
D.OtherMatrixVars.Rm = D.Rm.(doa_field);
D.OtherMatrixVars.Idx = D.Idx.(doa_field);



D.RadComp = D.RadVel;

D.RangeBearHead(:,2) = D.Bear.(doa_field);

% delete extra fields
fn = intersect(fieldnames(D),{'dBear','RmsBrgErr','BrgDiff','RomsBrg', ...
    'RadVel','eigValues','U','V','Flag', ...
    'Bear'});

D = rmfield(D,fn);

D.LonLat = rangeBear2LonLat(D.RangeBearHead,D.SiteOrigin);

D = compute_heading(D);


% unpack some of the struct fields

fn = {'RunTime','Pwr','LR'};

for i = 1:numel(fn)
    D.(fn{i}) = D.(fn{i}).(doa_field);
end


end

function R = fix_multicol_fields(R)
% FIX MULTICOL FIELDS - dole out fields with multi columns for cat in time


% Fix some fields prior  the temporal concat

% % these only for seasondes
% if size(R.RadVel,2) == 3
R.Params1 = R.Params(:,1);
R.Params2 = R.Params(:,2);
R.Params3 = R.Params(:,3);
% end

R = rmfield(R,'Params');

% R.SNR3 = R.SNR(:,3);
% R.SNR2 = R.SNR(:,2);
% R.SNR = R.SNR(:,1);

R.FileName = {R.FileName};


end


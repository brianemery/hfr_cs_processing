function run_cs_processing
% RUN CS PROCESSING - cross spectra to radials driver script
%
% ARCHITECTURE IDEAS
%
% - could follow range processing
% - arbitrary array geometry, fft length, etc
% - arbitrary doa method 
% - future use of things like different detectors, MAP algorithm, ...
% - tests using Anthony's data for validation, but also could build some
%   edge case tests using simiulations (eg 0-360 transition, etc)
%
%  
% NOTES
% will need to generalize CS format (CS.a1, CS.A78, eg?). 
%
% DEFAULT BEHAVIOR
% - Outputs radials every 10 min

% Copyright (C) 2017 Brian Emery
%
% Last update 26-May-2017 10:34:13


% TO DO
% - guess I need to specify a header.txt file to get the range increment
%   since it's not generally available anywhere else. For now assume that
%   distance to first range cell is same as the spacing (given in rkm)
% - add some user feedback "processing ... blah"
% - Document processing settings: Tau, FOL method and settings, site,
%   flist, APM used, what else ...
% 
% - This is a draft of what will be a cs_processing project release
%   (hfr_cs_proc?) 
% - Need test data set ... one set of CSQ and R to compare new with prev
%   maybe compare with a codar version of the radial ...
%
% - I have a 100cm/s cutoff, but sometimes get larger Vr - why?


% NOTE: start parpool or matlabpool (see line 42 of doa_on_cs.m)
%
% NOTE2: sci1 on vici with the below settings is the defacto test ... need
% to encode it ...


% CONFIGURE

site = 'sni1'; %'sci1'; % same as site file name in radial file

% FOL user parameters (see imageFOL.m help)
user_param = [40 100 5]; % [stdev maxcurr snrcutoff]

% Main Data Folder (Has CS, RDL subdirs and APMs in it)
wd = ['/home/emery/data/' site];
% wd = ['/projects/error_covariance/data/' site ];

% Define CS file directory
csd = [wd '/csq/'];

% Define radial data output directory 
rwd = [wd '/RDLs/'];

% apm file to use ... need a better way to organize these ...
apm_file = [wd '/MeasPattern_sni1_20130831.txt']; %'/MeasPattern_sci1_20130815.txt']; % 

% Determine time (min) of CS to average, set to zero if no averaging (eg if
% processing CSS). 
tau = 70; 
cs_type = 'CSQ';


% % case ambiguity ...
% site = upper(site);







% PRE - PROCESS

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





% PROCESSING LOOP

for i = 1:numel(rtimes)  
    
    % load cs files based on what is needed to average together
    CSL = cs_load(flist,ftimes,rtimes(i),tau,CSL);    
    
        
    % ship removal step ...
    % this needs to happen before computing the average
    CSL = cs_ship_rm(CSL);
        
    % compute the average of the CS for processing
    CS = cs_average(CSL); 

    
        
    % add standard CS info 
    [CS.freqs,CS.Vrad,CS.fb] = getVelocities(CS.Header);
    
    % do this now, use the data in the doa subfunction
    CS.SNR = get_SNR(CS);

    
    
    % FOL processing (rows are range cells) 
    [~, FOregi, ~] = imageFOLs(CS.antenna3Self.',iBragg,v_incr,user_param);
    
    % convert FOregi to peakIdx
    peakIdx = convert_fol_ix(FOregi);
    
    %  pcolor(10*log10(CS.antenna1Self.'))
    %  shading flat; colorbar; caxis([-160 -80]);
    %  hold on
    % h= plot(FOregi(:,2),FOregi(:,1),'c.');
    %
    % keyboard, close all
    
    % Run DOA processing
    % (music or mle, music metrics, music error, etc)
    [MU,ML] = doa_on_cs(CS,APM,peakIdx);
    
    
    % detection step - currently uses codar's music parameters
    ML = apply_test_result(ML); 
    MU = apply_test_result(MU);
    
    % explicitly setting this here
    MU.ProcessingSteps = {'imageFOLs.m','music.m','run_param_test.m'};
    ML.ProcessingSteps = {'imageFOLs.m','mle_ap.m'}; %,'run_param_test.m'};
    
        
    % clean up radial struct meta data
    Rmu = get_radial_meta(MU,APM,rtimes(i),rkm);
    Rml = get_radial_meta(ML,APM,rtimes(i),rkm);
    
   
    % write the radial struct to a mat file
    save([rwd Rmu.FileName],'Rmu','Rml')
    
end




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


CS = ReadCS(fname);

[~,Vrad,~] = getVelocities(CS.Header);

n = length(Vrad)/2;


% get indecies of Bragg lines
[~,iBragg(1)] = min( abs( Vrad(1:n)) ) ;
[~,iBragg(2)] = min( abs( Vrad(n:end)) ) ; 

% adjust last result
iBragg(2) = iBragg(2) + n - 1;

% velocity increment m/s
v_incr = mode(diff(Vrad))/100;

% get distance to first range cell
rkm = CS.Header.distToFirstRangeCell;

end

function peakIdx = convert_fol_ix(FOregi)
% CONVERT FOL OUTPUT 
%
% FORegi -- range/dop velocity indices of all points marked as FO
%
% peakIdx - cell array numel = # range cells, with row index of peaks
%

rng = unique(FOregi(:,1));

peakIdx = cell(length(rng),1);

for i = 1:length(rng)
   
    peakIdx{i} = FOregi(FOregi==rng(i),2)';
    
end


end

function R = get_radial_meta(R,APM,ftime,rkm)
% CONVERT DOA STRUCT TO HFRP RADIAL STRUCT

R.Type = APM.Type(1:4);

R.SiteName = APM.SiteName;

R.SiteOrigin = APM.SiteOrigin;

R.SiteCode = 1;

% R.FileName = fileparts_name(fname);

R.TimeStamp = ftime; % fnames_to_times(fname,['CSS_' R.SiteName '_'],'yy_mm_dd_HHMM');

R.FileName = ['RDLm_' R.SiteName '_' datestr(R.TimeStamp,'yyyy_mm_dd_HHMM') '.mat'];

R.RangeBearHead = [];

% Generate Range 
R.RangeBearHead(:,1) = R.RngIdx * rkm;

% convert to bearing to ccwE
% Note from APM.README.BEAR_Units.BEAR_Units =  'degCWN'
R.RangeBearHead(:,2) = cwN2ccwE(R.Bear);

R.RangeBearHead(:,3) = NaN; % need this still 

% Generate LonLat
R.LonLat = rangeBear2LonLat(R.RangeBearHead,R.SiteOrigin);


R.RadComp = R.RadVel;



% % Populate U, V
% [R.U, R.V, R.Error, R.Flag] = deal(NaN(size(R.RadComp)));
% 
% % Need to compute heading ...
% % * HERE *
% 
% % From loadRDLfile, unchecked really 
% R.U = R.RadComp .* cosd(R.RangeBearHead(:,3));
% R.V = R.RadComp .* sind(R.RangeBearHead(:,3));

% % Add CS file name
% R.CSFileName = CS.FileName;


% final cleanup - dBear not used on real data
R = rmfield(R,{'Bear','RadVel','RngIdx','dBear'});


end




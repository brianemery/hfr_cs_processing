function [MU,ML] = doa_on_cs(CS,APM,peakIdx,K) 
% DOA ON CS - run DOA processing on CS data 
% [MU,ML] = doa_on_cs(CS,APM,peakIdx,K)
%
% A more general version of radial_from_cs.m, uses imageFOL.m to get first
% order, and is called by run_cs_processing.m
%
% Make sure CS have not been converted to dBm
%
% INPUTS
% CS struct, APM struct, peakIdx
% ... K snapshots for music error
%
% OUTPUTS
% radial structs MU and ML based on MUSIC and MLE-AP
%
% FROM doa_on_range_cell.m
% 
% SEE ALSO
% radial_from_cs.m
%

% Copyright (C) 2017 Brian Emery
%
% version 10-Apr-2017 14:56:35
% derived from earlier versions of radial_from_cs.m

% Check for test
if strcmp('--t',CS), test_case, return, end

% CONFIGURE

if nargin < 4
    % Infer the number of data snapshots (about 3.5?)
    K = CS.Header.averagingTimeMin/((CS.Header.fftLength/CS.Header.SwRfreqHz)/60);
end

% make sure units are *NOT* dbm
if isfield(CS,'Units') && strcmp('dBm',CS.Units)
    CS = cs_dbm2volts(CS);
end


% GET THE DOAS
parfor i = 1:length(peakIdx) %  %PARFOR

    [MU(i),ML(i)] = doa_on_range_cell(CS,APM,peakIdx{i},i,K);
    
end


% Concat the data
MU = struct_cat(1,MU);
ML = struct_cat(1,ML);




end

function [MU,ML] = doa_on_range_cell(CS,APM,peakIdx,rdx,K)
% DOA ON RANGE CELL - custom local version
% [MU,ML] = doa_on_range_cell(CS,APM,peakIdx,rdx)
%
% Make sure CS have not been converted to dBm
%
% INPUTS
% ... K snapshots for music error
% ... peakIdx is just a vector here
%
% OUTPUTS
% radial structs MU and ML based on MUSIC and MLE-AP


% INIT OUTPUTS
% matrix of single and dual bearing solutions for both methods
% index refers to the index of the APM
ML = doa_struct(length(peakIdx)); 
MU = ML;

% put this here temporarily
MU.Type = 'music';
ML.Type = 'mle_ap';

% get array matrix for music error (SeaSondes)
A = get_array_matrix(APM);


% loop over peak indicies
for f = 1:length(peakIdx);
    
    fbin = peakIdx(f);
    
    % build covariance matrix     
    C = make_cov(CS,fbin,rdx);
    
 
    % RUN DF METHODS
    
    % Max liklihood - alt projection NEED TO CHECK CODE - rm ends?
    ix{1,1} = mle_ap(A,C,1);   
    ix{2,1} = mle_ap(A,C,2);   
    
    % Max liklihood
    % [ix,~] = mle(APM,C);
    
    % Music
    % keyboard
    % [~,sdx,ddx,tr,P] = music(APM,C);
    [~,mx,D,V] = music(A,C,APM.BEAR);
    
    % .. could use a clean up
    [P,tr] = run_param_test(D,V,A,mx{2});
    
    % stack indecies in a row vector
    % account for possibility of length(dualIdx) = 1
    mx = [mx{1} mx{2}']; %[sdx ddx(:)'];
    ix = [ix{1} ix{2}'];
    
    % music
    MU.Bear(f,1:length(mx))   = APM.BEAR(mx);
    MU.RadVel(f,1:length(mx)) = CS.Vrad(fbin);
 
    % mle
    ML.Bear(f,:)   = APM.BEAR(ix);
    ML.RadVel(f,1:length(ix)) = CS.Vrad(fbin);
    
    
    % Just use MUSIC Parameter test for now ...
    MU.Dual(f) = tr;
    MU.Params(f,:) = P;
    ML.Dual(f) = tr;
        
    % Guess these should take V,D as inputs ...
    MU.Err(f,1) = music_error(A,C,APM.BEAR,K,1,mx(1));
    MU.Err(f,2:length(mx)) = music_error(A,C,APM.BEAR,K,2,mx(2:end));
        
    % Get SNR data
    MU.SNR(f,:) = [CS.SNR.antenna3Self(fbin,rdx) ...
                CS.SNR.antenna13CrossSp(fbin,rdx) CS.SNR.antenna23CrossSp(fbin,rdx)];
    ML.SNR(f,:) = MU.SNR(f,:);
    
    clear ix

    %     % PLOTTING
    %
    %     % MUSIC dual and single
    %     hu    = plot(MU.Bear(f,2:3),MU.RadVel(f,2:3),'go');
    %     hu(2) = plot(MU.Bear(f,1),MU.RadVel(f,1),'g*');
    %
    %     % MLE dual and single
    %     hx    = plot(ML.Bear(f,2:3),ML.RadVel(f,2:3),'ro');
    %     hx(2) = plot(ML.Bear(f,1),ML.RadVel(f,1),'r*');
    %
    %
    %     if f == 1,
    %         legend([hu(:)' hx(:)'],'MUSIC dual','MUSIC sngl','MLE dual','MLE sngl')
    %
    %     end
    %
    %     pause
end


% need this range index to get range later
MU.RngIdx = rdx*ones(size(MU.Bear));
ML.RngIdx = rdx*ones(size(ML.Bear));




end



function test_case

load /projects/error_covariance/data/doa_on_cs_test.mat

[CS.freqs,CS.Vrad,CS.fb] = getVelocities(CS.Header);

    CS.SNR = get_SNR(CS);

tic
[MU,ML] = doa_on_cs(CS,APM,peakIdx);

toc
% With PARFOR
% Elapsed time is 27.734820 seconds.
%
% Witout PARFOR
% Elapsed time is 98.712372 seconds.
keyboard

end

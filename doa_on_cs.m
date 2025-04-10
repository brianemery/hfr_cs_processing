function S = doa_on_cs(CS,APM,peakIdx,K,CFG) %,n,dmth) 
% DOA ON CS - run DOA processing on CS data 
% S = doa_on_cs(CS,APM,peakIdx,K,CFG)
%
% A more general version of radial_from_cs.m, called by run_cs_processing.m
% and run_rng_processing.m 
%
% Make sure CS have not been converted to dBm
%
% INPUTS
% CS struct, APM struct, peakIdx as a cell
% ... K snapshots for music error
% peakIdx - A cell array (with numel = (# range cells)), each containing the
%           row index points within the First Order peaks.
% CFG should have:
% Nemit, dmth, use_parfor,mus_params ...
%
% OUTPUTS
% doa struct
%
% FROM doa_on_range_cell.m, which it now calls
%
% SEE ALSO
% radial_from_cs.m (old)

% Copyright (C) 2017 Brian Emery
%
% version 10-Apr-2017 14:56:35
% derived from earlier versions of radial_from_cs.m
%
% updates Nov 2018 to merge with doa_on_range_cell.m

% TO DO
%  
% appear to sometimes get dual= true when there is only one velocity
% solution ... see vr_averaging_investigation.m for evidence, but code here
% could probably test that there is a 2nd velocity and mark it single
% otherwise


% Check for test
if strcmp('--t',CS), test_case, return, end


% CONFIGURE, CHECK INPUTS

if nargin < 4
    % Infer the number of data snapshots (about 3.5?) ... less desirable
    K = CS.Header.averagingTimeMin/((CS.Header.fftLength/CS.Header.SwRfreqHz)/60);
    
    % disp(['Computed value of K from Header info: ' num2str(K) ])
    
end

if nargin < 5
    % also set emitters to search for 
    CFG.Nemit = 2;
    CFG.dmth = {'mu'}; 
    CFG.use_parfor = false; 
    CFG.mus_param = [10 5 8];
end
        
n = CFG.Nemit;


% make sure units are *NOT* dbm
if isfield(CS,'Units') && strcmp('dBm',CS.Units)
    CS = cs_dbm2volts(CS);
end
    
% Pre-pre allocate, more or less
[S(1:numel(peakIdx))] = deal(doa_struct(1,n,n,CFG.dmth));





% GET THE DOAS WITH CONDITIONAL PARFOR
% [~,r] = system('hostname'); % r = 'parfor disabled'

if ~(CFG.use_parfor) %|| strncmp(r,'yourcomputer',12) ... strncmp(r,'ekman',5) && 
       
    for rc = 1:numel(peakIdx)  
        
        S(rc) = doa_on_range_cell(CS,APM,peakIdx{rc},rc,K,CFG); %,n,dmth);
        
    end
    
elseif CFG.use_parfor 
        
    parfor rc = 1:numel(peakIdx) % %PARFOR over different range cells
        
        S(rc) = doa_on_range_cell(CS,APM,peakIdx{rc},rc,K,CFG); %n,dmth);
                
    end
        
end


% Concat all the range cells
S = struct_cat(1,S);

% Document
try, S.ProcessingSteps = {'doa_on_range_cell','doa_on_cs.m'}; catch, end


end


function test_case
% TEST CASE
%
% dev test data from cs_read.m

% Get cs test filename
fn = '/m_files/test_data/compute_apm_from_csq/CSQ_cop1_08_12_06_205124.cs';  

% probably not a valid APM for this CS but good enough for the functional
% test ...
load /m_files/test_data/doa_on_range_cell.mat APM


% load subset of range cells and see how they compare
CS = ReadCS(fn); 

[CS.freqs,CS.Vrad,CS.fb] = getVelocities(CS.Header);

% FOL processing (rows are range cells)
user_param = [40 100 5];

% ... from subfunction to run_cs_processing ...
[~,Vrad,~] = getVelocities(CS.Header);

[~,~,dv] = getDopplerVelocities(CS.Header);

n = length(Vrad)/2;

% get approx indecies of Bragg lines ... in general
[~,iBragg(1)] = min( abs( Vrad(1:n)) ) ;
[~,iBragg(2)] = min( abs( Vrad(n:end)) ) ; 

% adjust last result
iBragg(2) = iBragg(2) + n - 1;

% velocity increment m/s
v_incr = dv/100; %mode(diff(Vrad))/100;

% get distance to first range cell
rkm = CS.Header.distToFirstRangeCell;

% ...

[~, FOregi, ~] = imageFOLs(CS.antenna3Self.',iBragg,v_incr,user_param);

% convert FOregi to peakIdx
peakIdx = convert_fol_ix(FOregi);


% check plot
rc = 25;
plot(1:512,real(10*log10(CS.antenna12CrossSp(:,rc))),'-'), hold on
        
plot(peakIdx{rc},real(10*log10(CS.antenna12CrossSp(peakIdx{rc},rc))),'*'), hold on


S = doa_on_cs(CS,APM,peakIdx);
 
 
keyboard




% OLDER SIMPLER TEST (AND TIMES, PERHAPS)

% % .. a low SNR case from who knows where ...
% load /m_files/test_data/doa_on_range_cell.mat

[CS.freqs,CS.Vrad,CS.fb] = getVelocities(CS.Header);

CS.SNR = get_SNR(CS);


plot(CS.freqs,10*log10(CS.antenna3Self))
hold on
plot(CS.freqs(peakIdx),10*log10(CS.antenna3Self(peakIdx)),'ro')

% doa_on_cs requires cell input for peakIdx
px = peakIdx;
peakIdx = {px};
peakIdx{2} = px; % build up for parfor testing ... must match nmber of rc's



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

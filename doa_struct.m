function R = doa_struct(r,n,m,fn)
% DOA STRUCT - create empty struct for DOA solutions
% R = doa_struct(r,n,m,dmth)
%
% r = number of rows for each field
% n = number of emitters that will be searched for
% m = size of cov matrix
% dmth = cell of doa methods ({'MU','ML','WM','WF','SM'};)
%
% NOTES
%
% 1) that this creates large matricies containing the single, dual, three,
% ... up to n DOA solutions, so these are intended as an intermediate data
% structure that could be transitioned to a Radial struct after applying
% the hypothesis testing (aka signal detection) step. 
% 
% 2) These are created such that each radvel (in a given range cell) is a
% row, and so it should be possible to vertically concat these to create a
% Radial struct for a single time. And then temporally concat them after
% signal detection step. 
%
% see doa_on_range_cell, and doa_on_cs.m
%
% EXAMPLE
% MU = doa_struct(length(peakIdx));
%
% Notes
% dBear formerly used by compute_doa_errors.m, see
%    clean_up_and_recompute_errors.m, which calls
%    compute_doa_errors_wrt_roms.m


if nargin < 2, n = 2; end

if nargin < 3, m = n; end

if nargin < 4, fn = {'MU','ML'}; end

fn = upper(fn);
    
% Expand for all possible emitters
mx = (n*(n+1))/2;

R = RADIALstruct(1);

% Usually Ideal or Meas so putting the APM type here is appropriate
R.Type = '';

% Clear these out
[R.SiteName] = deal('');

% Add these for DOA processing
[R.RadVel,R.Bear,R.SNR,R.dBear,R.Err, ...
                        R.Params,R.Apprch] = deal(NaN(r,mx));

% Snr matrix
R.SNR = NaN(r,3);

% Guess at the size of these
R.eigValues = NaN(r,m);

% EigenVectors and Cov matrix as a cell array
R.OtherMatrixVars.eigVectors{r,1} = [];
R.OtherMatrixVars.cov{r,1} = [];

% dont need these here
R.OtherMatrixVars = rmfield(R.OtherMatrixVars,{'ERSC','ERTC','ESPC','MAXV','MINV','SPRC','VFLG'});

% Needed by CS processing
R.RngIdx = ones(r,1);

% Expand eventually for arbitrary arrays .. NaNs here to prevent errors
R.Dual = NaN(r,1);

% Track DOA compute time
R.RunTime = 0;

R(1).README.Bear = '(Struct) Estimated DOA for Method';
R(1).README.RadVel = 'Radial Velocity';
R(1).README.Params = 'MUSIC parameters';
R(1).README.SNR = 'From get_SNR.m';
R(1).README.dBear = 'DEPRECATED (Struct) Estimated (True) Brg Error (RMS)';
R(1).README.Err = 'SN89 stddev (deg)';
R(1).README.RunTime = '(Struct) DOA Method Total Compute Time (s)';
R(1).README.LR = 'Likelihood Ratio (glrt.M)';
R(1).README.Pwr = 'Signal Power (signal_power.m)';

R.README.RmsBrgErr = 'RMS bearing error (Struct with fields for DOA algo)';
R.README.BrgDiff = 'Bearing Difference (ROMS vs DOA method) (Struct with fields for DOA algo)';
R.README.RomsBrg = 'Bearing for nearest ROMS current (cwN)'; % see roms_to_radvel.m and m_idist.m
R.README.Rm = 'Model covariance from glrt.m';
R.README.Idx = 'struct of doa method detection APM indecies';
R.README.K = 'Estimated number of independent data snapshots';

R(1).RADIAL_struct_version = 'doa_struct';


% MAKE DOA METHOD SUBSTRUCTS

% if x, 
%     fn = {'MU','ML','WM','WF','SM'};
% else
%     fn = {'MU','ML'};
% end
% 

% create containers size r x mx
for i = 1:numel(fn)
    S.(fn{i}) = NaN(r,mx);
end

[R.Bear,R.dBear,R.RmsBrgErr,R.BrgDiff,R.RomsBrg,R.Pwr] = deal(S);


% create runtime count
for i = 1:numel(fn)
    S.(fn{i}) = 0;
end

[R.RunTime] = deal(S);


% create detection containers
for i = 1:numel(fn)
    S.(fn{i}) = NaN(r,n);
end

[R.LR] = deal(S);


% storage for model cov from glrt
for i = 1:numel(fn)
    S.(fn{i}) = cell(r,n);
end

[R.Rm,R.Idx] = deal(S);


% this is a single value for the whole CS file
R.K = NaN;



% OLDER - MAYBE USERFUL?

% % add detection matricies
% [S.MU,S.ML,S.WM,S.WF,S.SM] = deal(NaN(r,n));
% 
% [R.LR,R.GM,R.LR2] = deal(S);
% 
% % storage for model cov from glrt
% [R.Rm.MU,R.Rm.ML,R.Rm.WM,R.Rm.WF,R.Rm.SM] = deal(cell(r,n));


% 
% % Specify both DOA method and array type
% S.MU.Type = [APM.Type ' MUSIC']; 
% S.ML.Type = [APM.Type ' MLE-AP'];
% S.WM.Type = [APM.Type ' WMUSIC'];
% S.WF.Type = [APM.Type ' WSF'];
% S.SM.Type = [APM.Type ' SML'];


end
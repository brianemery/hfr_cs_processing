function R = doa_struct(n,m)
% DOA STRUCT - create empty struct for DOA solutions
% R = doa_struct(n,m)
%
% n = number of rows for each field
% m = max number of possible emitters
%
% NOTES
% 1) that this creates large matricies containing the single, dual, three,
% ... up to n DOA solutions, so these are intended as an intermediate data
% structure that could be transitioned to a Radial struct after applying
% the hypothesis testing (aka signal detection) step. 
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
% dBear is used by compute_doa_errors.m

if nargin < 2
    m = 2;
end
    
% Expand for all possible emitters
mx = (m*(m+1))/2;

R = RADIALstruct(1);

% Usually Ideal or Meas so putting the APM type here is appropriate
R.Type = '';

% Clear these out
[R.SiteName] = deal('');

% Add these for DOA processing
[R.RadVel,R.Bear,R.SNR,R.dBear,R.Err, ...
                        R.Params] = deal(NaN(n,mx));


% Guess at the size of these
R.eigValues = NaN(n,m+1);

% EigenVectors as a cell array?

% Needed by CS processing
R.RngIdx = ones(n,1);

% Expand eventually for arbitrary arrays
R.Dual = false(n,1);

% Track DOA compute time
R.RunTime = 0;

R(1).README.Bear = '(Struct) Estimated DOA for Method';
R(1).README.RadVel = 'Radial Velocity';
R(1).README.Params = 'MUSIC parameters';
R(1).README.SNR = 'From get_SNR.m';
R(1).README.dBear = '(Struct) Estimated (True) Brg Error (RMS)';
R(1).README.Err = 'SN89 stddev (deg)';
R(1).README.RunTime = '(Struct) DOA Method Total Compute Time (s)'; 

R(1).RADIAL_struct_version = 'doa_struct';


% Hack to make these into structures

[S.MU,S.ML,S.WM,S.WF,S.SM] = deal(NaN(n,mx));

[R.Bear,R.dBear] = deal(S);

[S.MU,S.ML,S.WM,S.WF,S.SM] = deal(0);

[R.RunTime] = deal(S);




end
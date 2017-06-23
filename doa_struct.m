function R = doa_struct(n)
% DOA STRUCT - create empty struct for DOA solutions
% R = doa_struct(n)
%
% n = number of rows for each field
%
% see doa_on_range_cell, and others
%
% EXAMPLE
% MU = doa_struct(length(peakIdx));
%
% NEW! TO DO base this on radial struct ...

% Notes
% dBear is used by compute_doa_errors.m


R = RADIALstruct(1);
% 
% R = radial_struct_ext(fname,A)

% Clear these out
[R.Type,R.SiteName] = deal('');

% Add these for DOA processing
[R.Bear,R.RadVel,R.Params,R.SNR,R.dBear,R.Err,R.RngIdx] = deal(NaN(n,3));

R.Dual = false(n,1);

R(1).README.Bear = 'Estimated DOA';
R(1).README.RadVel = 'Radial Velocity';
R(1).README.Params = 'MUSIC parameters';
R(1).README.SNR = 'From get_SNR.m';
R(1).README.dBear = 'Estimated Brg Error (RMS)';
R(1).README.Err = 'SN89 stddev (deg)';



end
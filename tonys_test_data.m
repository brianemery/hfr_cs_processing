function [C,APM] = tonys_test_data
% TONYS TEST DATA - data for testing DOA methods
% [C,APM] = tonys_test_data;
%
% TEST DATA EXAMPLE FROM 
% Properties of HF RADAR Compact Antenna Arrays and Their Effect 
% on the MUSIC Algorithm by dePaolo and Terril
%
% DUAL! BEARING EXAMPLE FROM GETDOA.M
% 
%
% EXAMPLE
% Makes a figure showing the DOA function vs bearing just like fig 9 in the
% De Paolo and Terrill Scripps report
%
% SEE ALSO
% music.m

% Copyright (C) 2016 Brian Emery

% Create the covariance matrix from de Paolo's example:
C=[ 0.2162          0.0303-0.0090i  0.3170-0.0063i; ...
    0.0303+0.0090i  0.0436         -0.0091+0.0213i; ...
    0.3170+0.0063i -0.0091-0.0213i  0.5416];

% Create the idealized pattern 
APM = make_ideal_pattern(225, 0:5:360);


return
% TESTING CODE BELOW - USEFUL FOR COPY AND PASTE

ix = mle_test(APM,C); % !!! THIS WORKS
APM.BEAR(ix)

keyboard

% MUSIC SOLUTION
[DOA,singleIdx,dualIdx] = music(APM,C);

% instead of using getDOAs, use the fact that the known signal input
% bearing are 205 and 330 deg (from the reference)
dualIdx = [find(A.BEAR == 205) find(A.BEAR == 330)]; %APM?


% FIGURES
figure

subplot(211)
plot_doa(APM,DOA(:,1),singleIdx,'Single Bearing DOA function')
axis([0 360 -4 12])


subplot(212)
plot_doa(APM,DOA(:,2),dualIdx,'Dual Bearing DOA function')
axis([0 360 -5 30])

end
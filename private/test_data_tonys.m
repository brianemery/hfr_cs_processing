function [C,APM] = test_data_tonys
% TONYS TEST ... EXAMPLE FROM 
% Properties of HF RADAR Compact Antenna Arrays and Their Effect 
% on the MUSIC Algorithm by dePaolo and Terril
%
% DUAL! BEARING EXAMPLE FROM GETDOA.M
% 
% Makes a figure showing the DOA function vs bearing just like fig 9 in the
% De Paolo and Terrill Scripps report
%
% NOTE that he made this from a simulation with currents input at 205 and
% 330 degrees which MUSIC gets kind of wrong
%
% Lots of code from music.m
% 
% (APM,DOA,singleIdx,dualIdx)


% Create the covariance matrix from de Paolo's example:
C=[ 0.2162          0.0303-0.0090i  0.3170-0.0063i; ...
    0.0303+0.0090i  0.0436         -0.0091+0.0213i; ...
    0.3170+0.0063i -0.0091-0.0213i  0.5416];

% Create the idealized pattern 
APM = make_ideal_pattern(225, 100:0.1:359);



end
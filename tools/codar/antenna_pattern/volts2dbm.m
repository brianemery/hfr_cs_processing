function [dBm,phs] = volts2dbm(v2)
% VOLTS SQUARED TO DBM - convert seasonde signal to dBm
% [dBm,deg] = volts2dbm(v2)
%
% Also outputs phase in radians if given two output arguments (and inputs
% are complex)
% 
% See File_CrossSpectra.pdf from the COS documentation
%
% See also cs_volts2dbm.m for converting units on a 
% whole cross spectra file.
%
% NOTE on units
% File_CrossSpectra.pdf suggests volts^2, which is correct even though 
% SpectraPlotterMap says volts. 'dBm' is a measure of power referenced to
% 1mW, suggesting that the CS files (which are conjugate products of the
% FFT's of the voltage time-series) contain power. Hence volts^2 
%
% See also mag2db, and db2mag for non HFR related applications (they use 
% 20*log10), magPhase2RealImag, abs, angle.


% Copyright (C) 2011 Brian M. Emery
% 14 Oct 2010 - Verified vs SpectraPlotterMap 
%  6 May 2011 - rechecked and renamed, added test case and
%               output of phase angle

%  if strcmp('--t',v2), test_case, return, end  %<-- comment out for speed

% From File_CrossSpectra.pdf:
% Note to convert Antenna1,2,3 to dBm use: 
% 10*log10(abs(voltagesquared)) + (40.0 ? 5.8) 
% The 40.0 is an adjustment to conversion loss in the receiver. 
% The 5.8 is an adjustment to processing gain in SeaSondeAcquisition.
%
% It seems the above is in error, this checks out with
% SpectraPlotterMap (see test below):
dBm = 10.*log10(abs(v2)) - 40 + 5.8;

% skip atand for speed if nargout is one
if nargout > 1,
    
    % Output phase in radians 
    % phs = atan(imag(v2)./real(v2));
    % ... on the interval (-pi to pi) ... required to reproduce complex #s
    phs = atan2(imag(v2),real(v2));
end

end

% --------------------------------------------------------
function test_case
% TEST CASE

% FROM SPECTRA PLOTTER
% Some numbers from CSQ_cop1_08_12_06_202548.cs, range 14, doppler 269
% A1, A2, A2 (shows volts??), 

v2 = [7.6137e-9 9.8172e-9 7.8471e-8];
dB = [-115.4 -114.3 -105.3];

run_check(dB,round(volts2dbm(v2).*10)./10)

% A13, A23, A12 (complex numbers)
v2 = [4.216e-9 - 2.408e-8i; 2.769e-8 + 1.846e-9i; 9.214e-10 - 8.596e-9i ];
dB = [-110.3; -109.8; -114.8 ];
ph = [-80.1; 3.8; -83.9]*pi/180;

keyboard

[dBm,phs] = volts2dbm(v2);

% This tests functionality, ... see dbm2volts for reversibility and sig
% digits
run_check(dB,round(dBm.*1e3)./1e3)
run_check(round(ph.*1e2)./1e2,round(phs.*1e2)./1e2)


% Here's a second check from SpectraPlotterMap (v11.1.2)
% (also shows volts instead of v^2)
v2 = [1.1906e-3 3.1938e-4 5.2024e-4];
dB = [-63.4 -69.2 -67.0];

run_check(dB,round(volts2dbm(v2).*10)./10)

keyboard
end

% --------------------------------------------------------
function run_check(db,v2)

if isequal(db,v2)
    disp('test passed')
else
    disp('test not passed'), keyboard
end


end
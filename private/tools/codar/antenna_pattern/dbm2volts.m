function v2 = dbm2volts(dBm,phs)
% DBM 2 VOLTS SQUARED - convert seasonde signal to in dBm to volts
% v2 = dbm2volts(dBm,phs)
%
% Inputs in dBm, with optional second input including the phase 
% in radians (on interval -pi to pi, such that a complex number 
% is produced:
%
% 
% 
% See File_CrossSpectra.pdf from the COS documentation, and also 
% https://en.wikipedia.org/wiki/Complex_logarithm
% "Problems with inverting the complex exponential function"
%
% See also cs_volts2dbm.m for converting units on a 
% whole cross spectra file.
%
% NOTE on units
% SpectraPlotterMap suggests that the units on the signal is volts, while
% File_CrossSpectra.pdf suggests volts^2. File_TimeSeries.pdf refers to
% volts as well as signal power ... (see note from Bill)
%
% See also mag2db, and db2mag for non HFR related applications (they use 
% 20*log10), magPhase2RealImag, abs, angle.


% Copyright (C) 2012 Brian M. Emery
%
% From volts2dbm.m

% optional test
if strcmp('--t',dBm), test_case, return, end

% get the magnitude back
v2 = 10.^( (dBm + 40 -5.8 )./10);
    
if nargin > 1
    
    % use the phase to compute the complex output
    v2 = v2.*cos(phs) + v2.*1i.*sin(phs);
    
end

% From File_CrossSpectra.pdf:
% Note to convert Antenna1,2,3 to dBm use: 
% 10*log10(abs(voltagesquared)) + (40.0 ? 5.8) 
% The 40.0 is an adjustment to conversion loss in the receiver. 
% The 5.8 is an adjustment to processing gain in SeaSondeAcquisition.
%
% It seems the above is in error, this checks out with
% SpectraPlotterMap (see test below):
% 
% %   Y = DB2MAG(YDB) computes Y such that YDB = 20*log10(Y).
% 
%   Copyright 1986-2007 The MathWorks, Inc.
%   $Revision $ $Date: 2009/10/16 06:11:22 $
% y = 10.^(ydb/20);


% dBm = 10*log10(abs(v2)) - 40 + 5.8;
% 
% % Output phase in degrees
% deg = atand(imag(v2)./real(v2));

end

% --------------------------------------------------------
function test_case
% TEST CASE

% FROM SPECTRA PLOTTER
% Some numbers from CSQ_cop1_08_12_06_202548.cs, range 14, doppler 269
% A1, A2, A2 (shows volts??), 

v2 = [7.6137e-9 9.8172e-9 7.8471e-8];
dB = volts2dbm(v2);

run_check(round(dbm2volts(dB)*1e14)/1e14,v2)


% Here's a second check from SpectraPlotterMap (v11.1.2)
% (also shows volts instead of v^2)
v2 = [1.1906e-3 3.1938e-4 5.2024e-4];
dB = [-63.4 -69.2 -67.0];

run_check(dB,round(volts2dbm(v2).*10)./10)


% CHECK FOR REVERSIBILITY 
% complex to db and back 
% ... from volts2dbm.m:
% A13, A23, A12 (complex numbers)
v2 = [4.216e-9 - 2.4081e-8i; 2.769e-8 + 1.846e-9i; 9.214e-10 - 8.596e-9i ];


[dBm,phs] = volts2dbm(v2);

v2_ = dbm2volts(dBm,phs);


% This tests reversibility ... to about 7 decimal places
run_check(round(v2*1e13)./1e13,round(v2_*1e13)./1e13)


% ... with zero phases?
v2 = [1.1906e-3 3.1938e-4 5.2024e-4];

[dBm,phs] = volts2dbm(v2);

v2_ = dbm2volts(dBm,phs);

run_check(round(v2*1e13)./1e13,round(v2_*1e13)./1e13)

keyboard

end

% --------------------------------------------------------
function run_check(db,v2)

if isequal(db,v2)
    disp('test passed')
else
    disp('test not passed'), 
end


end
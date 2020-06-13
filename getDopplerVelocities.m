function [freqs,Vrad,dv] = getDopplerVelocities(p)
% GET DOPPLER VELOCITIES compute bin frequency and radial doppler velocity (cm/s)
% [freqs,Vrad,dv] = getDopplerVelocities(Header)
%
% This version is different from the MUSIC.M version, in that the ocean
% current velocity is NOT output (for that use getVelocities.m). This outputs the
% velocity indicated by the Doppler shifted signal, which would be the
% velocity of the wave plus the ocean current, or radial velocity of a ship
% for example.
%
% REFERENCE: 
% Lipa, B., and D. Barrick. "Least-squares methods for the extraction of
% surface currents from CODAR crossed-loop data: Application at ARSLOE."
% Oceanic Engineering, IEEE Journal of 8, no. 4 (1983): 226-253.
%
% INPUTS:
% Header, the header structure from ReadCS.M, or cs_read.m. Defaults to a
% reasonable mid-range (13 MHz) HFR if given no input.
%
% OUTPUTS:
% freqs in Hz, 
% Vrad in cm s^-1
% dv, the doppler velocity bin width (cm s^-1)

% TO DO
% verify this works for 5, 24 and 42 Mhz. Check out the results vs SpecMap
% Finish the test case
% - maybe process a CSA file and then use the radial velocities to check
% this (something where there would be no averaging)

% Copyright(C) 2011 Brian Emery
%   version 1 June 2008
%             June 2014 outputs dv instead of dv/2

if nargin > 0 && strcmp('--t',p), test_case; return, end



% DEFINITIONS

% Define g
%
% http://en.wikipedia.org/wiki/Gravitational_acceleration
% "a conventional standard value of exactly 9.80665 m/s2"
g = 9.80665;

% speed of light in m/s
c = 299792458;  
 



% VERIFY INPUTS

% set defaults if no struct given
if nargin < 1
    p.sweepUp   = 0;
    p.fftLength = 512;
    p.SwRfreqHz = 2;
    p.freqMHz   = 13.49;
    p.SwBWkHz   = 100.7080;
end


% check for empties or zeros
if isempty(p.SwRfreqHz) || p.SwRfreqHz ==0 
    p.SwRfreqHz = 2;
end
if isempty(p.freqMHz) || p.freqMHz ==0
    p.freqMHz = 13.49;
end
if isempty(p.SwBWkHz) || p.SwBWkHz ==0
    p.SwBWkHz = 100.7080;
end



% CALCULATIONS
%
% Relevant Doppler equation:
%
% v = df*c / 2f
%
% where f is the transmit frequency

% Compute center tx freq (checks out with SpectraPlotterMap)
% see also File_CrossSpectra.pdf. Units: MHz
if p.sweepUp
    txfreq = (p.freqMHz + ((p.SwBWkHz/1000)*.5));
else 
    txfreq = (p.freqMHz - ((p.SwBWkHz/1000)*.5));
end



% This is the DopplerResolutionHzPerBin, also from Header.txt:
% 0.00390625                   ! 9 ALTERNATE Doppler Resolution in Hz
df = p.SwRfreqHz/p.fftLength;


% % Ok, we need 512 points total which doesn't split evenly around a zero. Codar
% % puts the one extra on the right side. 512 bins should be centered on the 
% % velocity?
% %
% % This puts 256 bins on either side of zero Doppler. I think this is 'more
% % correct' than SpectraPlotterMap
% nhalf = p.fftLength/2;
% 
% freqs = df:df:(nhalf*df);
% 
% freqs = [-freqs(end:-1:1) freqs];
%
% This is the SeaSonde way
nhalf = p.fftLength/2;
freqs = df .* [-(nhalf-1):0 1:nhalf]; 

% Compute velocity associated with Doppler shift, convert to cm/s
Vrad = (freqs*c/(2*txfreq) )*(100/1e6);


% Compute radar wavelength in meters
wvLen = (c/txfreq) / 1e6;

% Compute the velocity increment in cm/s
%
% Relevant formulas:
% k0 = wave number = 2*pi/wvLen
% df = 2 * k0 * dv  --> dv = df/(2*k0) 
%
% Should be the same as (-Vrad(1)+Vrad(end))/512
% the 100 converts m to cm
dv = (wvLen/2) * 100 * p.SwRfreqHz/p.fftLength;

freqs = freqs(:);

% keyboard
% 
% % Compare to Belinda's way:
% Nfft = p.fftLength;
% belVel = (-floor(Nfft/2) + (1:Nfft)) * dv;
% plot(Vrad,belVel,'-o')
% grid




end

function test_case
% TEST CASE
% 
% RUN BASIC TEST
%
% Feb 2015 TO DO:
% compare with ComputeDVel.m
%
% NOTES from calc_dv.m (Codar's code)
% Nfft = 512; 
% % Belinda's way
% belVel = (-floor(Nfft/2) + (1:Nfft)) * dv;
% plot(Vrad,belVel,'-o')
% % Use this if using fftshift?
% vel = (-floor(Nfft/2) + (0:Nfft-1)) * dv;
% hold on
% plot(Vrad,vel,'-ro')



% Define test file
fn = '/m_files/test_data/getDopplerVelocities/CSQ_cop1_08_12_06_200428.cs';

dat = cs_read(fn,':',1);
    
[freqs,Vrad,dv] = getDopplerVelocities(dat.Header);
     

% Get random sample from SpectraPlotterMap 
% These from CSQ_cop1_08_12_06_200428.cs 
% codar's index + 1 to convert to matlab. v is in m/s
i = [1 140 233 256 416  509];
v = [-11.11 -5.05 -1.00 0.00 6.97 11.02];

test_case_checks(Vrad(i)./100,v)



% Define test files
fn = '/m_files/test_data/getDopplerVelocities/CSQ_AGL1_12_05_04_172230.cs';

dat = cs_read(fn,':',1);
    
[freqs,Vrad,dv] = getDopplerVelocities(dat.Header);



% Get random sample from SpectraPlotterMap 
% These from CSQ_AGL1_13_03_07_001638.cs 
% codar's index + 1 to convert to matlab. v is in m/s
i = [1 51 280 511 782 999];
v = [-15.58 -14.06 -7.08 -0.03 8.23 14.85];

test_case_checks(Vrad(i)./100,v)






end

function test_case_checks(a,b)

rmsd = rmsdiff(a,b) * 100;

if rmsd < 3
    disp(['test ' inputname(2) ' ... ok (rmsdiff = ' num2str(rmsd) ', which is < 3 cm/s)'])
else
    disp(['test ' inputname(2) ' ... NOT ok'])
    keyboard
end


end

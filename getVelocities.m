function [freqs,Vrad,fb] = getVelocities(p)
% GET VELOCITIES.M compute bin frequency and radial ocean current velocity
% [freqs,Vrad,fb] = getVelocities(p)
%
% Referernces: Lipa and Barrick 1983 IEEE paper
%
% INPUTS
% p, the header structure from ReadCS.M, which contains:
%
%            dateTimeSec04: 3.3113e+09
%               kindOfData: 1
%                  siteStr: 'sci1'
%         averagingTimeMin: 4
%                  freqMHz: 13.4900         *
%                SwRfreqHz: 2               *
%                  SwBWkHz: 100.7080        *
%              sweepUpDown: 0               *
%              nRangeCells: 88
%           firstRangeCell: 0
%     distToFirstRangeCell: 1.4895
%                    dTHex: [4x1 double]
%         deleteRawSpectra: 0
%        overrideRemHeader: 0
%                fftLength: 512             *
%
% * denote variables used here!
%
% NOTE! USES REASONABLE DEFAULTS IF NO INPUT GIVEN!
%
% OUTPUTS
% freqs [Hz] (These from getDopplerVelocities.m) (COLUMN)
% Vrad  [cm s^-1]                                (COLUMN)
% fb    [Hz] Bragg Frequency (function of TX frequency)

% Copyright (C) 2011 Brian M. Emery
%       June, July 2008
%       June 2014 expand for other FFT lengths
% 


if nargin > 0 && strcmp('--t',p), test_case; return, end

% set defaults if no struct given
if nargin < 1
    p.sweepUp   = 0;
    p.fftLength = 512;
    p.SwRfreqHz = 2;
    p.freqMHz   = 13.49;
    p.SwBWkHz   = 100.7080;
end



% NEW NOTES - June 2014
% The following RDL file header key is wrong (SeaSonde v 6ish):
% DopplerResolutionHzPerBin: 0.001953125
% Values in Header.txt, TestReadCS.m and computed values compared with
% SpectraPlotterMap indicated the value is twice the above. It should be
% the sweep rate divided by the FFT length (eg 2/512 = 0.00390625);
%
% TO DO
% Devise a test (ship data?) to establish the accuracy here. 
% Determine if the center of the sweep is the tx freq to use, or if
%   something else would be better


% get the Doppler velocities
[freqs,~,~] = getDopplerVelocities(p);



% CODE FROM GET DOPPLER VELOCITIES

% Define g
%
% http://en.wikipedia.org/wiki/Gravitational_acceleration
% "a conventional standard value of exactly 9.80665 m/s2"
g = 9.80665;

% speed of light in m/s
c = 299792458;  


% Compute center tx freq (checks out with SpectraPlotterMap)
% see also File_CrossSpectra.pdf. Units: MHz
% ... now appears that down sweep has negative bandwitdh
%if p.sweepUp
    txfreq = (p.freqMHz + ((p.SwBWkHz/1000)*.5));
% else 
%     txfreq = (p.freqMHz - ((p.SwBWkHz/1000)*.5));
% end





% CODE FOR THIS MFILE

% Convert txfreq from MHz to Hz
txfreq = txfreq * 1e6;

% Compute radar wave number 
k0 = 2*pi*txfreq/c;

% Compute Bragg frequency (the frequency shift relative to txfreq due to
% the speed of the ocean wave). * SEP MFILE? *
fb = sqrt(2*g*k0)/(2*pi);

% Compute the ocean current speeds. This is done by useing the Doppler
% equation and removing the deltaF resulting from the ocean wave speed, in
% cm/s

n = p.fftLength/2;

% this is the frequency 'residual' after removing shift due to ocean wave speed
nfreqs = [freqs(1:n) + fb; freqs(n+1:end) - fb]; 

% This gets a difference answer than Spectral Plotter map, but it's close.
% UPDATE (Nov 2015): getDopplerVelocities now outputs adjusted freq, so
% this should be same as Spectra Plotter Map
Vrad = (nfreqs*c/(2*txfreq) )*(100);




% 
% 
% % compare with previous
% [pfreqs,pVrad] = previous_get_velocities(p);
% 
% keyboard
% 


return

end

function [freqs,Vrad] = previous_get_velocities(p)
% PREVIOUS VERSION - will rm later

% NOTES:
% This produces different values than CSPlot. The Radial Short data has no
% velocity values, as are shown in CSPlot. The radial short file data
% suggests that left and right sides are binned differently, which is what
% is done in this mfile, eg: my zero radial velocity indexies are 161 and
% 351 with values -1.67 and +1.67 cm/s - NOT zero. Details are illustrated
% in notesAboutVelocityBins.m. See also TestReadCS.m


% GET BRAGG FREQ
% Equations:
% wb = bragg Freq, k0= radar wave number, lambda = radar wave length
% wb = sqrt(2*g*k0); k0=2pi/lambda, lamda=c/f, -> k0=2pi*f/c
% --> fb = wb/2*pi = sqrt(2*g*ko)/2pi
% Here, I use the center TX freq. As given in File_CrossSpectra.pdf:
% "Sweep direction is up if non zero, else down.
% Center Tx freq is fStartFreqMHz+fBandwidthKHz/2 * -2^(bSweepUp==0)"
k0= 2*pi* (p.freqMHz-((p.SwBWkHz/1000)*.5)) *1000000/2.998e8; % cycles/m
braggFreqHz = sqrt(2*9.806*k0)/(2*pi);

binSpacingHz = p.SwRfreqHz/p.fftLength;



% ASSIGN DOPPLER FREQUENCY TO EACH BIN
% Re: 512: Seems that to get 512 points, codar eliminates -256. Also seems
% that codar defines the bragg bins as minimum radial velocity, then uses the
% freq bin spacing to get radial velocities from there. Force bins 161 and
% 351 to be at the bragg freq, and thus at radial velocity = minimum;
% 
% Other FFT length settings empirically determined
if p.fftLength == 512
    
    freqs=[([-160:0 1:95 ]*binSpacingHz) - braggFreqHz ...
           ([ -94:0 1:161]*binSpacingHz) + braggFreqHz ];
    
elseif p.fftLength == 256
    
    keyboard
    
elseif p.fftLength == 1024
    
    keyboard
    
else
    error(['codar:' mfilename ':fftLength not supported'])
end





% COMPUTE RADIAL VELOCITY FOR EACH BIN
% Equations:
% vrad = dw/2k0 (dw is delta freq from bragg)
% w = 2pi f -> vrad = 2pi * df / 2k0 = pi df/k0
% Force Left zero to -1.67, and right zero to +1.67 (in cm/s)
Vrad=[-1.67 + ( pi*(freqs(  1:end/2  )  + braggFreqHz)/k0 * 100 )...
       1.67 + ( pi*(freqs((end/2)+1:end)- braggFreqHz)/k0 * 100)]; % < -- why add then subtract braggFHz?

% Check numbers:
% "and a 50 cm/s current speed gives dw/2pi = 0.085 Hz" for 25.4Mhz
% Vrad=pi*[.085]/(2*pi*25.4e6/2.99e8)* 100
% = 50 cm/s
%
% delta freq = dw = 2k0*v(brg) ; v(brg) is radial velocity
% bragg freq in hz = wb/2pi ;
% eg at 25.4Mhz, wb = .514 Hz ...
%
% radar wavelength in meters: 2*pi/k0

Vrad = Vrad(:); freqs = freqs(:);
end

function test_case
% TEST CASE
% 
% RUN BASIC TEST



% Define test file
fn = '/m_files/test_data/getDopplerVelocities/CSQ_cop1_08_12_06_200428.cs';

dat = cs_read(fn,':',1);
    
[~,Vrad] = getVelocities(dat.Header);
     

% Get random sample from SpectraPlotterMap 
% These from CSQ_cop1_08_12_06_200428.cs 
% codar's index + 1 to convert to matlab. v is in m/s
i = [1 140 233 256 416  509];
v = [-693.8 -88.2 317.0 417.2 279.9 685.0];

test_case_checks(Vrad(i),v(:))



% Define test files
fn = '/m_files/test_data/getDopplerVelocities/CSQ_AGL1_12_05_04_172230.cs';

dat = cs_read(fn,':',1);
    
[~,Vrad] = getVelocities(dat.Header);



% Get random sample from SpectraPlotterMap 
% These from CSQ_AGL1_13_03_07_001638.cs 
% codar's index + 1 to convert to matlab. v is in m/s
i = [1 51 280 511 782 999];
v = [-860.2 -707.8 -9.4 695.1 125.3 787.1];

test_case_checks(Vrad(i),v(:))





end

function test_case_checks(a,b)

rmsd = rmsdiff(a,b);

if rmsd < 3
    disp(['test ' inputname(2) ' ... ok (rmsdiff = ' num2str(rmsd) ', which is < 3 cm/s)'])
else
    disp(['test ' inputname(2) ' ... NOT ok'])
    keyboard
end


end
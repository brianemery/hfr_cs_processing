function  snr = get_SNR(CS,idx)
% GET_SNR - signal to noise ratio for cross spectra using Codar method
% snr = get_SNR(CS,absIdx)
%
% Noise level is computed using cs_get_noise_level.m, an implementation of the
% Codar method. NOTE this previously used Self spectra to compute
% SNRs, now it uses this:
%
% fn = {'antenna3Self','antenna13CrossSp','antenna23CrossSp'};
%
% INPUTS
% CS        - data structure (from cs_read.m) in [volts^2] or dBm
% absIdx    - (Optional) the ABSOLUTE range-freq index of points to compute SNR on
%             
%
% OUTPUTS
% snr       - a matrix containing the SNR (col 1 is A33, 2 is A13, 3 is A23), 
%             and rows correspond to points in idx or the whole freq range
%
% NOTE ON OUTPUTS! If absIdx is not provided, this code will output a
% struct containing SNR for each antenna.  
%
% EXAMPLE
%
% 

% Copyright (C) 2009-2010 Brian Emery 
% 16 Sept 2009
% 15 Nov 2010 - updates for getNoiseLevel.m compatiblitity, matricized,
%               requires absolute index input, added subfunction
%  6 May 2011 - added test case, minor refactoring
%  9 Jun 2011 - added detection of input units
% 18 Aug 2014 - encorp'd cs_get_noise_level, refactored

% TO DO
% finalize verification of these SNR numbers with SpectraPlotterMap




% INITIALIZE AND INPUT CHECKS

% Optionally run test case
if strcmp(CS,'--t'), test_case, snr =[]; return, end

% Check for indexing
if nargin < 2, 
    
    idx = ':';
    
%     % convert matrix size to an absolute index
%     [r,c]=size(z); address = reshape(1:(r*c),r,c);
    
end



% Convert to dBm
CS = cs_volts2dbm(CS);

% Define fields of CS file
% fn = {'antenna1Self', 'antenna2Self', 'antenna3Self'};
fn = intersect(fieldnames(CS),{'antenna3Self','antenna13CrossSp', ...
                               'antenna23CrossSp','a33','a13','a23', ...
                               'a0303','a0103','a0203',});


% Compute noise levels for whole CS
CS = cs_get_noise_level(CS);


% COMPUTE SNR
%
% for whole CS

for i = 1:numel(fn)
    
    % compute snr after expanding noise scalars
    CS.SNR.(fn{i}) = CS.(fn{i}) - (ones(size(CS.(fn{1}),1),1)*CS.noiseLevel.(fn{i})) ;
    
end



% CREATE OUTPUTS

if ~ischar(idx)
    
    % preallocate
    snr = NaN(length(idx),numel(fn));
    
    for i = 1:numel(fn)
        
        % note the use of absolute indexing
        snr(:,i) = CS.SNR.(fn{i})(idx);
        
    end
    
else
    
    % I suspect this code won't be used much and that if I need the SNR for
    % the whole think I'd just output the SNR struct. So do that!
    
    
%     % preallocate
%     snr = NaN(length(idx),numel(fn));
%     
%     for i = 1:numel(fn)
        
        % note the use of absolute indexing
        snr = CS.SNR;
        
%     end
    
    
    
end


% -- CODE BELOW OLD METHOD FOR COMPARISON --

% GET SNR

% Use subfunction to get SNR for each point in idx
% for i = 1:numel(fn)
%     
%     
%     % Compute snr
%     snr(1:length(idx),i) = get_antenna_snr(CS.(fn{i}),freqs,idx);
%     
% end
%     


end

% --------------------------------------------------------
function test_case
% TEST CASE

% FROM SPECTRA PLOTTER
% Some numbers from CSQ_cop1_08_12_06_202548.cs, range 14, doppler 270
% (range, doppler start at zero) SNR's are:
% [27.3 32.7 35.5] (NF -142.3 -147.4 -140.9)

csqDataDir = '/m_files/test_data/cs_filter/good_cop_case/';

CS = cs_read([csqDataDir 'CSQ_cop1_08_12_06_202548.cs'],15);

snr = get_SNR(CS,270); 


% Uses numbers from previous run, checked vs spectraplottermap
run_check([37.3   33.4   36.2],round(snr.*10)./10)



% 2ND TEST WHOLE CS
CS = cs_read([csqDataDir 'CSQ_cop1_08_12_06_202548.cs']);

snr = get_SNR(CS); 

load /m_files/test_data/get_SNR/get_snr.mat SNR

% Uses numbers from previous run, checked vs spectraplottermap
run_check(SNR,snr)


return

% PLOT MAP OF SNR TO VERIFY
h = surf(SNR.antenna3Self');
view(2)
set(h,'EdgeColor','interp')
caxis([0 35])
colorbar

end

% --------------------------------------------------------
function run_check(db,v2)

if isequal(db,v2)
    disp('get_SNR.m: test ... ok')
else
    disp('get_SNR.m: test ...NOT ok'), keyboard
end


end




% OBSOLETE CODE ...

% This is here until cs_get_noise_level is incorp'd
function nBar = getNoiseLevel(css,freqs)
% GET NOISE LEVEL - compute noise level for each range cell
% nBar = getNoiseLevel(css,freqs,fb)
% 
% Compute Mean noise level using method detailed in
% FirstOrder_Settings.pdf. Output units: dBm
%
% INPUT
% css   - cross spectra matrix (eg CS.antenna1Self) in *dBm* - or convert
%         nBar to dDm later ...
% freqs - frequencies of velocity bins from getDopplerVelocity.m 
% fb    - (Optional) Bragg freq, needed for nFFT other than 512
%

% Copyright (C) 2009-2011 Brian M. Emery
% June 2008
%   Nov 2010 - reworked for matrix inputs
% 6 May 2011 - added test case, minor refactoring
% 9 Jun 2011 - lots of monkeying with the bins used to try to get it to
%              produce the same output as SpectraPlotterMap w/o success
% 16Feb 2012 - checked for .' vs ' effects






% Optionally run test case
if strcmp(css,'--t'), test_case, return, end

% this is totaly hardwired and wont work for fft lengths other than 512
if length(freqs) == 512
    
    % Columns are ranges, rows are frequencies
    % ideal bragg freqs are at 161 and 351 (row indicies)
    nse = [css(freqs >-0.96 & freqs <= freqs(161)-0.33,:); ...
           css(freqs < 0.96 & freqs >= freqs(351)+0.33,:) ];    
    
elseif length(freqs) == 1024
   
    disp('Warning: assuming BW = 25 hKz in getNoiseLevel.m')
    
   % generalize the above a bit 
    nse = [css(freqs >-0.4902 & freqs <= -fb - 0.1650,:); ...
           css(freqs < 0.4902 & freqs >=  fb + 0.1650,:) ];   
       
else
    disp('getNoiseLevel.m: need to support this FFT length')
    keyboard

    
end


   
% set Inf's to NaN
nse(isinf(nse)) = NaN;

% Get a standard deviation and mean of the noise region for each range cell
[nBar,stdev,~,~,~,~] = stats_noNaN(nse.');


% Recompute mean without points 3 stdev's out. Take the mean for each range
% cell, abs, add 3stdevs, then make it into a matrix with rows size nse
nse( abs(nse) > ones(size(nse,1),1)*(abs(nBar.')+(3*stdev.')) ) = NaN;

% Recompute
nBar = mean_noNaN(nse.'); nBar = nBar(:).';

end


% --------------------------------------------------------
function snr = get_antenna_snr(CS,freqs,ix)
% GET ANTENNA SNR 
% abstracted for one element of CS file

% Get noise level (dBm) for each range cell
n = getNoiseLevel(CS,freqs);

% Matrixize noise (such that ix applies to Noise also)
n = ones(size(freqs,1),1)*n;

% Compute SNR
snr = CS(ix) - n(ix);

keyboard

end


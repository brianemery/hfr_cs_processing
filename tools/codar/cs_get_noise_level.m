function  CS = cs_get_noise_level(CS)
% CS GET NOISE LEVEL - compute noise level for each range cell in a CS file
% CS = cs_get_noise_level(CS)
%
% *For backward compatibility, see subfunctions to eg get_SNR.m
% 
% Compute Mean noise level using method detailed in
% FirstOrder_Settings.pdf. Output units: dBm
%
% INPUTS
% CS        - data structure (from cs_read.m) in [volts^2] or dBm
%
% EXAMPLE
% % Define file, read it in, compute noise level
% fn = '/m_files/test_data/getDopplerVelocities/CSQ_AGL1_12_05_04_172230.cs';
%     
% CS = cs_read(fn); 
% 
% CS = cs_get_noise_level(CS);


% Copyright (C) 2014 Brian M. Emery
% 
% Version 15-Aug-2014 21:20:33

% ALGORITHM
% from qaqcCombinedDocument.pdf:
%"the mean ... [is] found for spectral points between about +/- 0.7 - 0.96
% Hz (10 frequency bins from the edge of the spectra to the ideal bragg
% frequencies +.33 Hz)".
%
% See also: Codar's FirstOrder_Settings.pdf



% INPUT CHECKS

% Optionally run test case
if strcmp(CS,'--t'), test_case, return, end


% Check for older useage (logical 3rd or 4th input)
if ~isstruct(CS)
    
    disp([mfilename ': input must be struct'])
    return
    
end

% Check and convert volts^2 to dBm
CS = cs_volts2dbm(CS);


% DEAL WITH MULTIELEMENT STRUCTS
if numel(CS) > 1
    
    [CS.freqs,CS.Vrad,] = deal([]);
    
    [CS.noiseLevel,CS.noiseStdev] = deal(struct([]));
   
   for i = 1:numel(CS)
       CS(i) = cs_get_noise_level(CS(i));
   end
   return
end



% INITIALIZE  

% Compute frequencies
[CS.freqs,CS.Vrad,fb] = getVelocities(CS.Header);


% Define fields of CS file ... old and new
fn = intersect(fieldnames(CS), ...
      {'antenna3Self','antenna13CrossSp','antenna23CrossSp', ...
      'antenna1Self', 'antenna2Self', 'antenna12CrossSp', ...
      'a33','a13','a23','a11','a22','a12', ...
      'a0303','a0103','a0203'});


% Document
CS.noiseLevel(1).README = 'Antenna noise level for each range cell (dBm)'; 
CS.noiseStdev(1).README = 'Stdev of noise for each range cell (dBm)'; 

% LOOP OVER CS FIELDS

% Use subfunction to get noise for each rangecell
for i = 1:numel(fn)
        
    % Compute noise level
    [CS.noiseLevel.(fn{i}),~,CS.noiseStdev.(fn{i})] = get_noise_level(CS.(fn{i}),CS.freqs,fb);
    
        
end
    


end

function [nBar,ix,stdev] = get_noise_level(css,freqs,fb)
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



% GENERALIZED FOR ALL FFT LENGTHS
ix = find( freqs > (0.96 * min(freqs)) & freqs <= (-fb - 0.33*max(freqs)) | ...
           freqs < (0.96 * max(freqs)) & freqs >= ( fb + 0.33*max(freqs)) );


nse = css(ix,:);
      
   
% set Inf's to NaN
nse(isinf(nse)) = NaN;

% Get a standard deviation and mean of the noise region for each range cell
[nBar,stdev] = stats_noNaN(nse.');


% Recompute mean without points 3 stdev's out. Take the mean for each range
% cell, abs, add 3stdevs, then make it into a matrix with rows size nse
nse( abs(nse) > ones(size(nse,1),1)*(abs(nBar.')+(3*stdev.')) ) = NaN;


% Recompute
nBar = mean_noNaN(nse.'); nBar = nBar(:).';




end

function test_case
% TEST CASE 
%
% Tests:
% cell list/multi element struct
% noise floors computed are same as previous (Aug 2014) implementation

fname = {'/m_files/test_data/cs_filter/good_cop_case/CSQ_cop1_08_12_06_202548.cs'; ...
    '/m_files/test_data/getDopplerVelocities/CSQ_AGL1_12_05_04_172230.cs'; ... % 1024 with 25 khz
    '/m_files/test_data/getDopplerVelocities/CSQ_ptc1_13_05_27_000615.cs';}; %256


CS = cs_read(fname); 

CS = cs_get_noise_level(CS);



% test vs previously compute data
load /m_files/test_data/cs_get_noise_level/cs_get_noise_level.mat NS

fn = {'antenna3Self','antenna13CrossSp','antenna23CrossSp'};



for j = 1:numel(NS)
    for i = 1:numel(fn)
        rslt(i,j) = isequal_tol(CS(j).noiseLevel.(fn{i}),NS(j).(fn{i}));        
    end    
end

run_check(rslt)

 keyboard


% TEST NEW IMPLEMENTATION
% 
% 
% % TEST NFFT > 512
% % Define test file
% 
% dat = cs_read(fn,':',1);
%     
% [freqs,Vrad,fb] = getVelocities(dat.Header);
% 
% nBar = getNoiseLevel(volts2dbm(dat.antenna1Self),freqs,fb);
% 
% keyboard
% 
% volts2dbm(nBar)
% 
% keyboard
% 
% 
% keyboard
% 
% 



return


% GENERATE CHECK PLOTS 
% These commented out but can provide a useful test of how things are
% working


% Define range cell
rc = randi(size(CS.antenna1Self,2),3); 

test_case_fig(rc(1),CS)
test_case_fig(rc(2),CS)
test_case_fig(rc(3),CS)

keyboard



% MAKE TEST DATA

for i = 1:3
    NS(i) = CS(i).noiseLevel;
end

for i = 1:3
    NS(i).FileName = CS(i).FileName;
    
end


save /m_files/test_data/cs_get_noise_level/cs_get_noise_level.mat NS 


end

function test_case_fig(rc,CS)

figure

% plot spectra
H = cs_plot(CS,rc); hold on


% add the noise levels
fn = fieldnames(CS.noiseLevel);
fn = setdiff(fn,{'README'});

for i = 1:3
    h = plot(H(i).ax,[-1 1],CS.noiseLevel.(fn{i})([1 1]),'-m');
    set(h,'LineWidth',2)

end



end

function run_check(db)

if all(db)
    disp([mfilename ': test ... ok '])
else
    disp([mfilename ':test ... NOT ok']), keyboard
end


end

function tf = isequal_tol(a,b)

% should be a funciton ...
tol = 1e3;

tf = isequal(round(tol*a)/tol, round(tol*b)/tol);

end

function [xbar,stdev] = stats_noNaN(x)
%
% ** LOCAL VERSION to speed calcs **
%
% From: 
% STATS_NONAN
% [xbar,stdev,hi,lo,median,n] = stats_noNaN(x)
% Computes the mean, std, and ranges of each row of a
% matrix after removing NaN's. Also returns number of
% non-NaN values used in the calc.
%
% See also mean_noNaN.m

% Copyright (C) 1999-2010 Brian M. Emery
% Brian Emery 13Sept99
% try-catch added 3apr09

% % check for test case
% if strcmp('--t',x), test_case, return, end
% 


r =size(x,1);

% % if x is row vector, orient it so that it computes the
% % mean of the vector
% if r==1 || c==1
%    x=x(:)';
%    r=size(x,1);
% end

% create empty outputs
% [xbar,stdev,hi,lo,med,n]=deal(NaN(r,1));
[xbar,stdev] = deal(NaN(r,1));

% If there is a way to do this without a loop, I'd like to know what it is.
for i=1:r
        row = x(i,~isnan(x(i,:)));
        
       % n(i) = length(row);
    xbar(i) = mean(row);
   stdev(i) = std(row);
     % med(i) = median(row);

%   % when NaN's are put in, these return empty matricies and an error
%   try hi(i) = max(row); 
%   catch
%       hi(i)=NaN; 
%   end
%   
%   try lo(i) = min(row); 
%   catch
%       lo(i)=NaN; 
%   end

end



end
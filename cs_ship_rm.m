function CS = cs_ship_rm(CS)
% CS SHIP RM - remove ship signal from CS data
%
% INPUT
% CS  - multielement struct of Cross Sectra (CSQ) data ... uses a field
%       name 'ProcessingSteps' to determine if this has been run on the 
%       individual CS data files already
%
% REFERENCE 
% US Patent 5,361,072 Gated FMCW DF Radar and Signal Processing
% ... Barrick et al. 1994

% Copyright (C) 2017 Brian Emery
%
% Version 26-May-2017 13:37:12

% TO DO
% -  need to do everything non-db as the conversion to dBm loses the imag
%    component, or I need to keep that ...
% - Range dependent SNR or normalize ship signal relative to max to get ships
%   at long range ... need some cases of this to identify

% NOTE
% if a ship gets into the Bragg region, it will cause a small increase in
% the power, but not one that is so big that it will be detected. Also, I
% dont think that subtracting power from that bin would 'reveal' the signal
% from the currents. So there is likely to still be errors resulting from
% the ships going through Bragg. Or does it?  I could look at some data 
% (SCI) and find this. ... Actually, I think the residual might just be the
% short term signal due to ships ...
% One way around this might be to output a SNR above BG for each bin in
% first order region and then use that the identify suspicious signal
% 

% IDEA
% Could use the SNR of the BG to identify anomalous signal levels due to
% transient things like ships, and create CS-sized matrix of flags that would get
% propagated to the radial level

% check for test case
if strcmp('--t',CS), test_case, return, end

% check for previous application of this step - true if done before 
% ... need a better way to track this 
tf = strcmp('cs_ship_rm',[CS.ProcessingSteps]);

% use this to get index of CS structs to loop over
ix = 1:numel(CS); ix = ix(~tf);   

% set default number of CS to compute background 
% *NOTE* ... this code wont work for
% other values of n at the moment
n = 3;

% All this needs to happen in dBm ... at least for the ship stuff
CdB = cs_volts2dbm(CS);

% % preallocate for noise calc
% [CdB.noiseLevel,CdB.noiseStdev] = deal(struct([]));



for i = ix  %3:numel(CS)-1 % this is a problem for n~=3
    
    % get index of n nearest CS to use in average
    cx = get_index_of_n(1:numel(CS),i,n);
   
    

    
    % compute the 'background' eg average n CS together out of the available
    CSA = cs_average(CdB(cx)); % try divide by median also (e.g. Barnum, 1986 IEEE)
    % CSA = cs_median(CdB(cx));
    
    % subtract out average to get the background 
    BG = cs_subtract(CdB(i),CSA);
    
    % smooth in doppler space prior to SNR calculation to 
    % further separate noise from ships
    BS = cs_spatial_smooth(BG,5);
    
    
    % get snr to facilitate getting high values for removal
    snr = get_SNR(BS); % snr = snr.antenna3Self;
    
      
    
 
    %  subtract the short-term (residual) signal power 
    CS(i)  = cs_subtract_ship(CdB(i),BG,CS(i),snr);
    
    
    % document ... 
    CS(i).ProcessingSteps(1) = {'cs_ship_rm'};
    
end



end

function  m = get_index_of_n(m,i,n)
% for m integers, find the n nearest to i


% sort by the distance from i
[~,ix] = sort( abs(m-i) ); 

% reorder the m integers
m = m(ix);

% get the first n, excluding i
m = m(2:n+1);


end


function CS = cs_subtract_ship(CdB,BG,CS,SNR)
% CS SUBTRACT SHIP
% CS(i)  = cs_subtract_ship(CdB(i),BG,CS,snr);
% .. from the dB ...
%  just where the ship signal is found (using SNR cutoff) ... then
%  convert back to volts2 and insert into CS


% Variable cutoff - start with 10 dB - codar seems to use 12
cut = 10;

% Find Doppler bins with high residual SNR
tf =   ( SNR.antenna3Self >= cut & ...
        (SNR.antenna13CrossSp >= cut | SNR.antenna23CrossSp >= cut ) );


% subtract background from CdB, whole thing
CdB = cs_subtract(CdB,BG);


% now convert CdB back to volts, and update only relevant data in CS
CdB = cs_dbm2volts(CdB);

fn = cs_fieldnames(CS);


for i = 1:numel(fn)
    
    CS.(fn{i})(tf) = CdB.(fn{i})(tf);
        
end

end

function BG = cs_subtract(CS,CSA)
% All inputs in dB ...

BG = CS;

fn = cs_fieldnames(CS);

for i = 1:numel(fn)
    
    BG.(fn{i}) = CS.(fn{i}) - CSA.(fn{i});
    
    % try this for the phases - think this what I want
    % ... note that zeros pass throuth
    dp = CS.Phases.(fn{i}) - CSA.Phases.(fn{i});
    BG.Phases.(fn{i}) = atan2(sin(dp),cos(dp));
    
end



end

function AV = cs_median(CS)
% CS MEDIAN - compute median of fields in multi element CS structs
% CSA = cs_median(CS)
%
% Uses all the data in the multi-element input ... robust to empties
% 
% typically applied to CS after conversion to dB - more meaninful that way

% Copyright (C) 2017 Brian Emery


% get field names
fn = cs_fieldnames(CS);

% create output
AV = CS(1);

    
for i = 1:numel(fn)
    
    % need this to get the size
    X = cat(3,CS.(fn{i}));
    
    % compute median
    AV.(fn{i}) = median(X,3);
    
end



end

function BG = cs_spatial_smooth(BG,bins)
%
% quick and dirty and probably fine ... side lobes? 
% 
% figure,
% plot(1:512,BG.antenna3Self(:,15))



fn = cs_fieldnames(BG);

for i = 1:numel(fn)
    
    BG.(fn{i}) = mvave3(BG.(fn{i}).',bins).';
    
    % % apply to phases too? ... no, this only used for SNR which doesnt need
    % % phase info
    % BG.Phases.(fn{i}) = mvave3(BG.Phases.(fn{i}).',bins).';
    
end

% hold on
% plot(1:512,BG.antenna3Self(:,15))






end


function test_case
% TEST CASE
% 
% test case directory: /m_files/test_data/
%
% Get data together for this
% wd = '/projects/soo_apm/data/Data_COP1/for_figure/';
% 
% flist = get_file_list(wd,'CSQ*');
% 
% for i = 1:numel(flist)
%     CS(i) = ReadCS(flist{i});
% end
% 
% save /m_files/test_data/cs_ship_rm.mat CS
%
% NOTE that CS(5) below appears to show a ship going through the Bragg
% region, and this code not really affecting the signal there ... which
% means that this code mostly only detects ships that would cause large
% anomolous velocities that look like first order 

% % NOTES
% tic, CSA = cs_average(CS); toc
% Elapsed time is 0.013143 seconds.


% Get data
load /m_files/test_data/cs_ship_rm.mat CS

[CS(1:numel(CS)).ProcessingSteps] = deal({''});


% before
H = cs_plot_map(CS(5));  % 4, 6 good to look at also
title(H.axes(1),'BEFORE')
  

% remove ships
CS = cs_ship_rm(CS);


% after
H = cs_plot_map(CS(5));
title(H.axes(1),'AFTER')



keyboard


% older code for debugging



s = cs_plot_map(BG);
for x = 1:3, caxis(s.axes(x),[-40 40]), end

title('BACKGROUND')

keyboard

H = cs_plot(CS(i),15)
figure, H = cs_plot(BG,15)
figure, H = cs_plot(CSA,15)

keyboard


% EXTRA
% numbers for checking things ... span a big range of SNR and are complex

fn = '/m_files/test_data/compute_apm_from_csq/CSQ_cop1_08_12_06_205124.cs';
CS = cs_read(fn);
x  = CS.antenna23CrossSp(359:366,10)



end



% delete these eventually

function BG = cs_zero_out(BG,snr,cut)



% strict-er snr criteria
tf = ~( snr.antenna3Self > cut & ...
       (snr.antenna13CrossSp > cut | snr.antenna23CrossSp > cut ) );
    
 
fn = cs_fieldnames(CS);


for i = 1:numel(fn)
    
    BG.(fn{i})(tf) = 0;
    
end

% 
% figure, H = cs_plot_map(BG);
% for x = 1:3, caxis(H.axes(x),[-30 30]), end
% 
% keyboard

end


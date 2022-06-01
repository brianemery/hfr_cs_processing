function [peakIdx,Alims,absIdx,abstf] = getFirstOrder(dat,Vrad,user_param)
% GET FIRST ORDER.M - Get the indecies of the first order region
% [peakIdx,Alims,absIdx,nBar]=getFirstOrder(CS,Vrad,user_param)
%
% For an entire Cross Spectra data file. see FirstOrder_Settings.pdf,
% or SpectraPlotter Map
%
% INPUTS
% CS    Cross Spectra Struct from cs_read.m or ReadCS.m
% Vrad  Ocean current velocities from getVelocities.m 
% user_param These are [Flim Fdown NoiseFac Nsmooth Currmax] (see FirstOrder_Settings.pdf)


% 
% OUTPUTS
% peakIdx - cell array numel = # range cells, with row index of peaks
% Alims   -  can be compared with Alims.txt, output by
%           SpectraPlotterMap.m, and also found in the Processing folder.
% absIdx - absolute indexing of peakIdx
% abstf  - logical array of CS matrix regions defined by peakIdx (for
%           plotting)
%
% Data used from Typical header.txt file:
% 100 4          !11 Maximum Current Limit, #Pts Smoothing For FirstOrderCalc
% 15.0 1 100.0   !12 Factor Down Peak limit 1st order Radials ,0 = no 2nd order, Factor down 1rst order Waves
% 7.5 4.0        !15 Radials:Factor down peak nulls, Noise factor,
%
% (lower noise accept data closer to the noise level
%
% NOTES:
% Header => CSPlot translations:
% 'signal above noise' = Above noise factor (eg 6db =>4)
% ' peak drop off' = Below peak factor
% ' allowed peak null' = above null factor
%
% antilog base b of y = b^y ... so V^2 is 10^(dbm/10)
%
% OUTLINE:
% 1) elimiate noisey data using est of noise in the wings and a threshold
%       threshold = noise level * noise fact
% 2) smooth the spectrum (monopole self spectrum)
% to the left and  right halves of the spectrum:
% 3) find the nulls using the smoothed spectrum, defining the first order region
%       Search for nulls on the perifery, at amax/fdown
% 4) elimiate points that are too far below peak energy, ie if 
%       elimiate if P(v^2)<amax/flim
% 5) eliminate doppler freqs greater than curr max
%
% SEE ALSO simulateDivergenceProblem.m


% 2014 IDEAS
% Use logical indexing, turn this into one big logical statement, output a
% logical firstOrder matrix, same size as CS data, 1 if first order. Seems 
% like this could be greatly simplified - but then it always does! Convert
% everything to same size matrix: velocity, smoothed CS, factor stuff ...
%
% "Could probably re-implement this as a logical/find type operation"
% TO DO
% - this is a total CF. Should just be one big logical statement
% - do everything in dB from the outset
% - apply to all range cells (test for CS formatting)
% - Then: screen by max current, snr, 
% - use my SNR code
%
% see subfunction to bearing_uncertainties.m



% Copyright (C) 2008-2010 Brian Emery
% Version 1, June 2008
% Version 2 
%   fixed: mv ave centered
%          zero crossing better
%          reproduced v1 results
% Version 2.1 
%   checked for ' vs .' effects (untested)


% 
%
%
% OLDER
% TO DO
% - Output row and absolute index versions of both alims and stuff in peak
% - Add examples using each (see test_case)
% - better code to seperate left and right sideds of bragg peak
% - encorp use of volts2dbm.m (test after)
%
% TO TRY
% - smooth over only 4 points?
% - experiment with positioning of step one. ( 1) elimiate noisy data?)
%    2) smooth ENTIRE? spectrum, 3)find nulls using limited search, 
% - Alims different than final peak in far range cells?




% Optionally, run test case (set plts = 2 below)
if strcmp(dat,'--t')
    test_case, return
end

% --------------------------------------------------------------
% USER SETTINGS (eg from header.txt) 
% --------------------------------------------------------------

if nargin < 3
    % DEFAULT SETTINGS
    % flim defines how far down from the spectra peak to accept
    header.flim        = 15.0; %15.85;   % these are in Volts^2 (they get converted to dBm) (Default 15.0)
    header.fdown       = 7.5; %5.01; %    % these are in Volts^2, defines where to start looking for nulls
    header.noiseFactor = 4.0; %10.0; %    % " " " " in Volts^2, defines how far above noise to accept
    header.currmax     = 100;    % cm/s
    
    % number of points to smooth for finding nulls (looks like what is actually
    % used is this number *2+1! (see SpectraPlotterMap)
    header.nsm = 4;
    
else
    
    header.flim        = user_param(1);
    header.fdown       = user_param(2);
    header.noiseFactor = user_param(3);
    header.currmax     = user_param(5);
    header.nsm         = user_param(4);

    
end



% % CUSTOM FOR PWS ANALYSIS
% disp([mfilename ':FIRST ORDER SETTINGS NON-STANDARD'])
% header.flim        = 15.85;    % 15.0;   % these are in Volts^2 (they get converted to dBm)
% header.currmax     = 150;      % 100 % cm/s
% header.fdown       = 5.01;     % 7.5;   % these are in Volts^2
% header.noiseFactor = 10;       % 4.0; % " " " " in Volts^2
% % number of points to smooth for finding nulls (looks like what is actually
% % used is this number *2+1! (see SpectraPlotterMap)
% header.nsm = 2; % 4

% Plotting verbosity switch: 0 for none,1 for one last plot, 2 for all of them
plts = 0; % 0, 1 or 2

% Range cells to plot:
RC2Plot = 1:53;% plots all of them %:45 50:52];%30:3:55];
%  disp('RCs 43-45 have 3 1 and 4 data points!')
%  disp('RCs 50-52 should have NO data!!'), disp('show which side of the peak stuff comes from ...')


% --------------------------------------------------------------
%  INITIALIZE
% --------------------------------------------------------------

% get monopole self spectrum, in dbm ... and loops 1, 2 self 
dat = cs_volts2dbm(dat);

% assign short names
mss= dat.antenna3Self; %10*log10(abs(dat.antenna3Self));% +40-5.8; % see File_CrossSpectra.pdf
A1 = dat.antenna1Self; %10*log10(abs(dat.antenna1Self));
A2 = dat.antenna2Self; %10*log10(abs(dat.antenna2Self));

% initialize variables
bins=size(mss,1);
Left.Idx=1:round(bins/2)-10;
Right.Idx=round(bins/2)+10:bins;
aboveNoiseIdx=cell(size(mss,2),1);
peakIdx=cell(size(mss,2),1);
Alims(1:size(mss,2),1:4)=NaN;
absIdx =[];

% use subfunction to get smoothed spectrum which is really step 2, but is
% best done all at once.
MssSmooth = mvAve(mss.',(header.nsm*2)+1,1).'; 


% Compute Noise levels (dBm)
dat = cs_get_noise_level(dat);

nBar = dat.noiseLevel.antenna3Self;


% --------------------------------------------------------------
% LOOP OVER EACH RANGE CELL
% --------------------------------------------------------------

for i = 1:size(mss,2)
    
%  % "Steps 1-2 are applied to the entire monopole self - spectrum"

    % ------------------------------
    % STEP 1.
    % loop over the range cells to find the points above the noise
    % "the noise level in the wings of the spectrum is calculated adn data
    % eliminated below this noise level times a factor 'noisefact'. For example
    % if noiseFact = 4, data will be accepted that is heigher than 4x the noise
    % floor"
    % Here I'm using the method for computing noise level that is based on the
    % method of computing SNR detailed in qaqcCombinedDocument.pdf. Need a ref.
    %
    % ... Some investigation suggests that this in applied to three of the
    % spectra ...
    %     nBar=getNoiseLevel(mss(:,i),freqs);
    %     nBarA1=getNoiseLevel(A1(:,i),freqs);
    %     nBarA2=getNoiseLevel(A2(:,i),freqs);
    %     % get points above noise
    %     n = find(mss(:,i)> (nBar+10*log10(header.noiseFactor)) & ...
    %         A1(:,i)> (nBar+10*log10(header.noiseFactor)) & ...
    %         A2(:,i)> (nBar+10*log10(header.noiseFactor)) );
    
    % Conversation with D. Barrick suggests that it should be done this
    % way: (self or cross spectra?)
    n = find(mss(:,i) > (nBar(i) + 10*log10(header.noiseFactor)) & ...
              A1(:,i) > (nBar(i) + 10*log10(header.noiseFactor))   & ...
              A2(:,i) > (nBar(i) + 10*log10(header.noiseFactor)) );
    aboveNoiseIdx{i}=n;

    checkPlotStep1(mss(:,i),i,RC2Plot,n,nBar(i),header,dat,plts)
    
    
    
    % ------------------------------
    %  STEP 2: 
    % Smooth the spectrum 
    % Note that in CSPlot it appears that the moving average is 
    % header.nsm*2+1 points ... hmmm
    mssSmth = MssSmooth(:,i);

    checkPlotStep2(mss(:,i),i,RC2Plot,n,nBar(i),mssSmth,dat,plts)
    

    
    % ------------------------------    
    % "steps 3-5 are applied seperately for the +ve and -ve Doppler halves of
    % the monopole self spectrum

    
    
    % ------------------------------
    %  STEP 3:
    % Find Nulls
    % "find the point with the maximum power (amax) [in dbm?], and the
    % surrounding points amax/fdown" this is where the search for the nulls
    % begins ...
    % Runs subfunction on halfs of each range cell.
    
    % LEFT SIDE: (idx is relative to leftIdx, thus make it absolute)
    Left = findNulls(mssSmth(Left.Idx),Vrad(Left.Idx),Left,header,i); 
    % RIGHT SIDE
    Right = findNulls(mssSmth(Right.Idx),Vrad(Right.Idx),Right,header,i);
    
    checkPlotStep3(mss,mssSmth,i,RC2Plot,Left,Right,dat,plts)
  
    
    
    % ------------------------------
    %  STEP 4:
    % Points are eliminated if power in V^2 is less than amax/flim
    % NOTE: apply to unsmoothed data?
    % LEFT SIDE:
    Left=limitSpectralRange(mss(:,i),Left,header);
    % RIGHT SIDE
    Right=limitSpectralRange(mss(:,i),Right,header);
    
    checkPlotStep4(mss,mssSmth,i,RC2Plot,Left,Right,dat,plts)



    % ------------------------------
    % STEP 5:
    % Use currmax to define frequency window to elimiate large false
    % currents. Vrad,Half,header
    % LEFT SIDE:
    Left=limitUsingCurrMax(Vrad,Left,header);
    % RIGHT SIDE
    Right=limitUsingCurrMax(Vrad,Right,header);
    
    checkPlotStep5(mss,mssSmth,i,RC2Plot,Left,Right,dat,plts)

    % keyboard
    
    % ------------------------------
    % MAKE ALIMS
    % keep the indicies and make Alims:
    % here I want to keep points that have passed all above peak finding
    % tests, along with the final 'aboveNoise' test in step 1.
    Left.PeakIdx=intersect(Left.PeakIdx,aboveNoiseIdx{i});
    Right.PeakIdx=intersect(Right.PeakIdx,aboveNoiseIdx{i});
    peakIdx{i}=[Left.PeakIdx(:)' Right.PeakIdx(:)']; 
    
    % create the Alims
    if ~isempty(Left.PeakIdx)
        Alims(i,1:2)=[min(Left.PeakIdx)  max(Left.PeakIdx)];
    end
    if ~isempty(Right.PeakIdx)
        Alims(i,3:4)= [min(Right.PeakIdx) max(Right.PeakIdx)];
    end
    
    % Final check plot, post application of Step 1 ...
    checkPlotStep6(mss,i,RC2Plot,peakIdx,Alims,dat,plts)
    
    
    % COMPUTE ABSOLUTE INDEX
    % Convert index of each range cell into an absolute index
    % Method using just points in the peak
    % if ~isempty(peakIdx{i})
    %    absIdx = [absIdx sub2ind(size(mss),peakIdx{i},i*ones(size(peakIdx{i})))];
    % end
    %
    % METHOD USING ALIMS
    % Use the Alims to get *all inclusive* 1st order region (might be
    % bigger than just the points above the noise)
    idx = [Alims(i,1):Alims(i,2) Alims(i,3):Alims(i,4) ];
    idx = idx(~isnan(idx));
    %
    % % This previous method appears to be wrong ... 
    %     if ~isempty(idx)
    %         absIdx = [absIdx sub2ind(size(mss),idx,i*ones(size(idx)))];
    %     end
    
end % end range cell loop




% CREATE OUTPUTS

% logical array
abstf = false(size(A1));

% absolute index of peakIdx
absIdx = [];

% loop over range cells again
for i = 1:size(mss,2)

    abstf(peakIdx{i},i) = true;
    
    absIdx = [absIdx sub2ind(size(mss),peakIdx{i},i*ones(size(peakIdx{i})))];
    
end


end



 
% Step Subfunctions
function out = mvAve(dataIn,window,cutoff)
% MVAVE4.M
% flt=mvave(data,i,j);
% Filter data using a centered window moving average. Where:
% 	'data' is the data to be filtered (row vector or matrix)
% 	'i' is a number which gives the size of the filter
% 	'flt' is the filtered data
% 	'j' is the required amount of data within the window
% 		(j must be a whole number, j<i)
%
% This function computes a moving average for each row if 'data'
% is a matrix. Uses technology developed for mean_noNaN.m to
% compute moving average without NaN's, such that a mean is computed
% from j or more amount of data within the moving average window.
%
% The output data is the same size as the inputs, with the first
% i/2 and last i/2 consisting of the unfiltered data.

% Version 1 Written by Brian Emery, sometime in 1995, using Libe
% Washburn's idea.
% Version 4 - 5oct1999 by Brian Emery
% Version here customized to do even or odd moving average

% get fractional cutoff
fract=cutoff./window;

% get size of data2, make space for outputs
[rows,cols]=size(dataIn);
out=dataIn;

for k=1:rows

    % name the current row 'data'
    data=dataIn(k,:);

    % define the matrix to be averaged
    D(1:window,length(data)+(window-1))=NaN;

    % use the loop to build the matrix to be averaged
    for m=1:window
        D(m,m:length(data)+(m-1))=data;
    end

    % filter the data using mean_noNaN to take the average
    flt=mean_noNaN_sub(D,fract);

    % keep only the part you want
    flt=flt(window:end-(window-1)).';

    % basically padding with actal data ...
    out(k,(1:length(flt))+floor(window/2))=flt;

%         % optional check plot
%         h1=plot(1:length(data),data,'-b.'); hold on,
%         h2=plot(1:length(out(k,:)),out(k,:),'ro');
%         title('MvAve.m Subfunction Check PLot')
%         legend([h1 h2],'data In',[num2str(window) 'pt Filtered'])
%         keyboard

end
end
 
function Half = findNulls(rcellSmth,Vrad,Half,header,i)
% Run this on each half of one range cell (of the monopole)
% separately.
%
% Step 3. Find nulls between 1st and 2nd order
% "find the point with the maximum power (amax) [in dbm?], and the
% surrounding points amax/fdown" this is where the search for the nulls
% begins
%
% Outputs the indecies of the NULLS (ie ALIMS?).

% Brian Emery Giugno, 2008



% NOTES
% Lipa and Barrick 83 has a version of peak picking which the currently
% used method is probably a derivative of. They take the diff of the
% smoothed spectrum and find the highest and lowest values which hopefully
% occur at the gap between 1st and 2nd order.


% USER SETTINGS ----------------------------
% set tolerance for finding zero crossing
tol=0.2; 
% set debug plot switch
plts=false;

% ------------------------------------------


% Find AMAX. Look only in the expected bragg region.
expectedIdx = find(Vrad(:) < header.currmax & Vrad(:) > -header.currmax);
[amax,amaxIdxRel]=max(rcellSmth(expectedIdx)); %amax=max(rcell); %
amaxIdx=expectedIdx(amaxIdxRel); 


% Using SPectraPlotterMap, this appears to be the part for "we dont want to search
% within the first order region", ie, the peak drop off setting
% DEFINE THRESHOLD amax/fdown
thresh=amax-10*log10(header.fdown);

% TAKE DIFF AND FIND POINTS BELOW THRESHOLD (looking for nulls)
diffRcell=diff(rcellSmth);

% only consider points below the threshold. NOTE: i'm eliminating points
% near the edges using curr max. Codar does not appear to do this, but I
% Think I have to do it because of my moving average filter. Maybe points
% near the ends should be filled with NaN instead?
idx=find(rcellSmth<thresh & Vrad(:) < 2*header.currmax & Vrad(:) > -2*header.currmax);
% idx=find(rcell<thresh & Vrad(:) < 2*header.currmax & Vrad(:) > -2*header.currmax);
idxLeft=idx(idx<amaxIdx); 
idxRight=idx(idx>amaxIdx); % these indicies apply to rcellSmth

if ~isempty(idxLeft) && ~isempty(idxRight) %<-- should be an or? or run as a subfxn on both separet
    % Get Last and first Zero crossings
    % First, get max and min indicies which define the left and right
    % halves of the peak (not to be confused with the left and right halves
    % of the spectrum)

    % find the zero crossing of the left half
    zerosLeft=find(diffRcell(idxLeft)>-tol & diffRcell(idxLeft)<tol);
    
    % find the zero crossing of the right half
    zerosRight=find(diffRcell(idxRight)>-tol & diffRcell(idxRight)<tol);
     
    % error catch if these are empty. Run iterations multiplying the
    % tolerance (this max 5 before, now set to max 10). This just means
    % that the closest diff is a large distance away. For the CSQ's this
    % might just mean the data is noisier?
    it=2;
    while isempty(zerosLeft) && it < 30
        if plts, disp(['zerosLeft empty! RC=' num2str(i) ' fixing ...']), end
        tol2=tol*it;
        zerosLeft=find(diffRcell(idxLeft)>-tol2 & diffRcell(idxLeft)<tol2);
        it=it+1; % if it>30, disp('it>30'),  end %keyboard,
    end
    it=2;
    while isempty(zerosRight) && it < 30
        if plts, disp(['zerosRight empty! RC=' num2str(i) ' fixing ...']), end
        tol2=tol*it;
        zerosRight=find(diffRcell(idxRight)>-tol2 & diffRcell(idxRight)<tol2);
        it=it+1; if it>30,  disp('it>30'), end %keyboard, end
    end
    
    try
        % Then get abs indecies of smoothed: nulls and points in peak
        nullIdx=[idxLeft(zerosLeft(end)) idxRight(zerosRight(1))];
        
    catch
        nullIdx = [];
        
    end
        
    % debug plots?
    if plts && i==13
        figure % freq space ...
        h1=plot(1:length(rcellSmth),rcellSmth,'-b.'); hold on
        h2=plot(idxLeft,rcellSmth(idxLeft),'r.');
        h3=plot(idxRight,rcellSmth(idxRight),'m.');
        h4=plot(nullIdx,rcellSmth(nullIdx),'gs');
        title(['Range Cell ' num2str(i) ': Find Null Debug: Freq Space'])
        legend([h1 h2 h3 h4],'smth spec','below thresh (L)','below thresh (R)','NULLs')
        keyboard

        figure % diff
        h1=plot(1:length(diffRcell),diffRcell,'-b.'); hold on
        h2=plot(idxLeft,diffRcell(idxLeft),'r.');
        h3=plot(idxRight,diffRcell(idxRight),'m.');
        h4=plot(idxLeft(zerosLeft),diffRcell(idxLeft(zerosLeft)),'go');
        h5=plot(idxRight(zerosRight),diffRcell(idxRight(zerosRight)),'go');
        title('Find Null Debug: Diff ')
        legend([h1 h2 h3 h4 h5],'diff spec','below thresh (L)','below thresh (R)', ...
                                'near Zero(L)','near Zero(R)')        
        pause
    end

    % put output into structure. Use absolute indicies
    Half.NullIdx=Half.Idx(nullIdx);
    
else % case when nulls cant be found? Such as when too little signal   
    % put output into structure. Use absolute indicies
    Half.NullIdx=[];
    
end

% output more data
Half.DiffIdx=Half.Idx(idx);
Half.Thresh=thresh;
Half.Amax=amax;
Half.diffRcell=diffRcell;

end
 
function Half = limitSpectralRange(rcell,Half,header)
% Step 4. Limiting the spectral range (applied after finding nulls)
% points are eliminated if power in V^2 is less than amax/flim
% NOTE: apply to unsmoothed data?
% Applies to one half of one range cell

% Define Outputs
Half.PeakIdx=[];

if ~isempty(Half.NullIdx)
    % Get indicies from null positions for this peak
    indx=Half.NullIdx(1):Half.NullIdx(2);

    % Apply limits.
    pdx=find(rcell(indx) >= Half.Amax-10*log10(header.flim));

    if ~isempty(pdx)
        Half.PeakIdx=indx(pdx);
    end
end
end
 
function Half = limitUsingCurrMax(Vrad,Half,header)
% Step 5. use currmax to define frequency window to elimiate large false
% currents
if ~isempty(Half.PeakIdx)
    vdx=find(Vrad(Half.PeakIdx)<=header.currmax & Vrad(Half.PeakIdx)>-header.currmax);

    if ~isempty(vdx)
        Half.PeakIdx=Half.PeakIdx(vdx)';
    else
        Half.PeakIdx=[];
    end
end
end

 
function abar = mean_noNaN_sub(ad,f)
% MEAN_NONAN
% abar=mean_noNaN(ad,f)
% Computes the mean of each COLUMN(!!) of a matrix after removing NaN's.
% This subfunction is a modified version which outputs the mean
% only when there are more than a changeable amount of data points
% on which to base the mean.

% Based on Code stolen from Erik Fields
% Updated 5Oct1999 Brian Emery

% Find the NaN's, set them to zero for the sum part of the mean
% calc, then set the zeros in the denominator to NaN. It's ugly,
% just figure it out...
%ad=a';
i=isnan(ad);
ad(i)=0;
sm=sum(~i); j=find(sm==0); sm(j)=NaN;
abar=(sum(ad)./sm).';
%if sum(~i)==0
%   abar=NaN;
%else
%   abar=(sum(ad)./sum(~i))';
%end

% abar is mean of all data without NaNs. Make rows of abar which
% have too many NaN's into NaN's. Fraction is the percent data in
% each row. If this fraction is less than that givn by
% f, the mean is set to NaN.
[r,c]=size(i);
fraction=sum(~i)./r;
abar(find(fraction<f))=NaN;
abar(find(isnan(fraction)))=NaN;

% force output to be a column vector
abar=abar(:);

end



 
% Check plotting Subfunctions
function checkPlotStep1(mss,i,RC2Plot,n,nBar,header,dat,plts)
% Step 1 Check Plots


if plts >1 % && ismember(i,RC2Plot)
    
    disp(['Running getFirstOrder.m Check Plotting, RC=' num2str(i)])
    
    figure
    h1=plot(1:length(mss),mss,'-b.'); hold on
    h2=plot([0 512],[nBar nBar],'k'); 
    h3=plot(n,mss(n),'ro');
    [h4,h5]=checkPlotAddLoopSpectra(dat,i);
    title(['STEP 1: Points above noise+factor, Range Cell ' num2str(i)])
    legend([h1 h2 h3 h4 h5],'spectrum', 'noise level', ...
        ['points above noise+factor = ' num2str(nBar+10*log10(header.noiseFactor))], ...
        'Loop1Self','Loop2Self')
     keyboard %pause
end
end
 
function checkPlotStep2(mss,i,RC2Plot,n,nBar,mssSmth,dat,plts)
% Step 2 check plots. builds on figure from Step 1
if plts >1 && ismember(i,RC2Plot)
    figure
    h1=plot(1:length(mss),mss,'-b.'); hold on
    h2=plot([0 512],[nBar nBar],'k');
    h3=plot(n,mss(n),'ro');
    h4=plot(1:length(mssSmth),mssSmth,'-m.');
    title(['STEP 2: ... with smoothed spectrum, Range Cell ' num2str(i)])
    legend([h1 h2 h3 h4],'spectrum', 'noise level', ...
        'points above noise+factor ','smoothed spectrum')
    pause
end
end
 
function checkPlotStep3(mss,mssSmth,i,RC2Plot,Left,Right,dat,plts)
%  'mss' here needs to be indexed for the range
% cell. Also:
%    Left.NullIdx=LeftIdx(idx);
%    Left.DiffIdx=LeftIdx(diffIdx);
%    Left.Thresh=thresh;  
if plts >1 && ismember(i,RC2Plot)
    figure
    h1=plot(1:length(mss(:,i)),mss(:,i),'-b.'); hold on
    h2=plot(1:length(mssSmth),mssSmth,'-m.'); keyboard
    h4=plot([Left.DiffIdx Right.DiffIdx],mssSmth([Left.DiffIdx Right.DiffIdx]),'g.');
    h3=plot([Left.NullIdx Right.NullIdx],mssSmth([Left.NullIdx Right.NullIdx]),'ro');
    a=axis; da=(a(2)-a(1))/2;
    h5=plot([a(1) a(2)-da],[Left.Thresh Left.Thresh],'k:');
    h6=plot([a(1)+da a(2)],[Right.Thresh Right.Thresh],'k--');
    title(['STEP 3: Left and Right with Nulls, RC=' num2str(i)])
    
    legend([h1 h2 h3 h4 h5 h6],'spectrum','smoothed spectrum','Nulls',...
        'diff points', 'Left amax/fdown','Right amax/fdown')    
    pause
    
    % Check the differencing is working
    figure
    h1=plot([Left.Idx(1:end-1) Right.Idx(1:end-1)],...
        [Left.diffRcell(:)' Right.diffRcell(:)'],'-b.'); hold on    
    h2=plot(Left.DiffIdx,Left.diffRcell(Left.DiffIdx),'m.');  
    h3=plot(Right.DiffIdx-1,Right.diffRcell(Right.DiffIdx-Right.Idx(1)-1),'r.');  
    legend([h1 h2 h3],'diff all points','pts below thresh (L)','pts below thresh (R)')
    pause
end
end
 
function checkPlotStep4(mss,mssSmth,i,RC2Plot,Left,Right,dat,plts)
% Check plot for step 4. show how spectra was limited
if plts >1 && ismember(i,RC2Plot)
    figure
    h1=plot(1:length(mss(:,i)),mss(:,i),'-b.'); hold on
    h2=plot(1:length(mssSmth),mssSmth,'-m.');
    h4=plot([Left.DiffIdx Right.DiffIdx],mssSmth([Left.DiffIdx Right.DiffIdx]),'g.');
    h3=plot([Left.NullIdx Right.NullIdx],mssSmth([Left.NullIdx Right.NullIdx]),'ro');
    h5=plot([Left.PeakIdx Right.PeakIdx],mss([Left.PeakIdx Right.PeakIdx],i),'c.');
    title(['STEP 4: Showing How Spectra Was Limited, RC=' num2str(i)])
    
    legend([h1 h2 h3 h4 h5],'spectrum','smoothed spectrum','Nulls',...
        'diff points', 'Peak after Amax/Flim')    
    pause
end
end
 
function checkPlotStep5(mss,mssSmth,i,RC2Plot,Left,Right,dat,plts)
% step 5: Show effect of CurrMax application
% Check plot for step 4. show how spectra was limited
if plts >1 && ismember(i,RC2Plot)
    figure
    h1=plot(1:length(mss(:,i)),mss(:,i),'-b.'); hold on
    h2=plot(1:length(mssSmth),mssSmth,'-m.');
    h4=plot([Left.DiffIdx Right.DiffIdx],mssSmth([Left.DiffIdx Right.DiffIdx]),'g.');
    h3=plot([Left.NullIdx Right.NullIdx],mssSmth([Left.NullIdx Right.NullIdx]),'r*');
    h5=plot([Left.PeakIdx Right.PeakIdx],mss([Left.PeakIdx Right.PeakIdx],i),'ks');
    title(['STEP 5: After Curr Max ... , RC=' num2str(i)])
    
    legend([h1 h2 h3 h4 h5],'spectrum','smoothed spectrum','Nulls',...
        'diff points', 'Peak after Amax/Flim, CurrMax')    
    pause
end

end
 
function checkPlotStep6(mss,i,RC2Plot,peakIdx,Alims,dat,plts)
% step 6: FInall plot, and data 'above noise' after step 1.
if plts >1 && ismember(i,RC2Plot)
    figure
    h1=plot(1:length(mss(:,i)),mss(:,i),'-b.'); hold on
    h2=plot(peakIdx{i},mss(peakIdx{i},i),'ks'); 
    alim=Alims(i,~isnan(Alims(i,:)));
    h3=plot(alim,mss(alim,i),'r*'); 
    title(['STEP 6?: After "above noise" Final Peaks!, RC=' num2str(i)])   
    legend([h1 h2 h3],'spectrum','Final Peak','Alims')    
    pause, close all
end

end
 
function [h1,h2]=checkPlotAddLoopSpectra(dat,i)
% add other antennas:

A1=10*log10(abs(dat.antenna1Self));
A2=10*log10(abs(dat.antenna2Self));
h1=plot(1:size(A1,1),A1(:,i),'.r-');
h2=plot(1:size(A2,1),A2(:,i),'.g-');

end

 
function test_case
% TEST CASE
% internal consistency check and a rudimentary check of Alims compared with
% data from SpectraPlotterMap

% test case to build: range cell with NO peaks, make sure this gets it
% right and downstream code can handle empties (March 2013)


% csq file to look at
csqNm = '/m_files/test_data/getFirstOrder/CSQ_Rfg1_10_07_21_060515.cs';

% Read it in ...
CS = ReadCS(csqNm,1);

% Get frequencies
[CS.freqs,CS.Vrad] = getVelocities(CS.Header);


% Run the test
[peakIdx,Alims,absIdx,abstf] = getFirstOrder(CS,CS.Vrad);



% This test currently not passing. Needs work to get it exactly same as
% codars. Super confusing
if isequal([156 179 330 348],Alims(20,:))
    disp('test SpectraPlotterMap consistency ... passed')
else
    disp('test SpectraPlotterMap consistency ... NOT passed')
end



% VISUAL CHECK 

% get stuff for plotting
[CS.freqs,CS.Vrad,dv] = getDopplerVelocities(CS.Header);

% Map plot
H = cs_plot_map(CS);

hx = H.axes;


test_case_plot(hx(3),CS,absIdx,peakIdx,Alims)


% PARTICULAR RANGE CELL
% This one seems to be difficult to get the whole peak due to it being
% split

CS = cs_volts2dbm(CS);

ix = peakIdx{10};

figure

H = cs_plot(CS,10);

hold on

plot(H(1).ax,CS.freqs(ix),CS.antenna13CrossSp(ix,10),'-mo')

plot(H(1).ax,CS.freqs(Alims(10,:)),CS.antenna13CrossSp(Alims(10,:),10),'b*')

keyboard

end


function test_case_plot(hx,CS,absIdx,peakIdx,Alims)
% TEST CASE PLOT
%
% add to each loop on the CS MAP

% get number of range cells
rc = size(CS.antenna1Self,2);


% make range and freq matricies for plotting
[R,F] = meshgrid(1:rc,CS.freqs);
%plot(F,R,'r.')


hold on

% Plot the absolute index of stuff in peak
plot(hx, F(absIdx),R(absIdx),'m.')

disp('checking absIdx (cyan) ... hit key to continue')


% plot the range cell index (these should be the same as absIdx)
for i = 1:size(Alims,1)
    plot(hx,CS.freqs(peakIdx{i}),i*ones(size(peakIdx{i})),'g.')
end

disp('checking range cell index (green) ... hit key to continue')



% plot the Alims
absAlims = [];

% convert Alims to absolute. Need to do this by hand b/c of the NaNs
for i = 1:size(Alims,1)

    rw = Alims(i,~isnan(Alims(i,:)));
    
    absAlims = [ absAlims; sub2ind(size(CS.antenna1Self), rw, i*ones(size(rw)))'];
    
end

plot(hx,F(absAlims),R(absAlims),'m.')

disp('checking Alims (magenta) ... hit key to continue')











end







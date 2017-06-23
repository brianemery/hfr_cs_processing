function R = radial_from_cs(flist,A) 
% RADIAL FROM CS.M - Get Radial data from Cross Spectra Using COS SeaSonde Method
% R = radial_from_cs(css_flist,APM)
%
% Reverse Engineer Codar's implementation of of MUltiple SIgnal
% Classification (MUSIC), to get radial files from CrossSpectra.
%
% SEE ALSO
% doa_on_cs.m, doa_on_range_cell.m, cs_processing
%
% NOTE
% this is sort of obsolete, see the above instead

% Copyright (C) 2014 Brian M. Emery
%       3 Jun 2008 Original implementation
%      21 Jul 2011 minor updates, modifications
%      17 Jun 2014 major refactor 
%
% Version 11-Jul-2014 11:11:17 (CSE = pi)


% TO DO
% move plots and tests to subfunction
% maybe extract music part and call this compute_radials or?
% keep info on single/dual bearing
% compute heading
% put plotting code somewhere (project toolbox?)
% 



% check for test case
if strcmp('--t',flist), test_case, return, end
 
% check ins, 
if ischar(flist)
    flist = cellstr(flist);
end



% Initialize RADIAL structure, one for each CS file
R = get_radial_struct(flist,A);


% Loop over each Cross Spectra File
for i = 1:numel(flist)
    
    R(i) = get_radial_from_cs(flist{i},R(i),A);
    
    
end
    



end

function R = get_radial_from_cs(fname,R,A)
% RADIAL FROM CS - get radial vel data from one CS file
%
% INPUTS 
% fname - CS file name (char), one element
% R     - Radial data struct (to be filled)
% A     - APM data struct


% Add other info to the RADIAL struct based on the CS file
NM = cosFileNameParts(fname);
R.FileName =['RDLm_' R.SiteName '_' NM.wholeDateStr '.ruv'];
R.TimeStamp = NM.TimeStamp;



% Load CS file
CS = ReadCS(fname);

% Compute radial velocities
[freqs,Vrad] = getVelocities(CS.Header);

% Get First order region
% peakIdx is a cell array with left and right indecies for each range cell.
[peakIdx,~,absIdx,~] = getFirstOrder(CS,Vrad);


% Expand data matricies based on CS data info
[R.RadComp, R.DOA, R.TestResults, ...
             R.MusicParamValues, R.freqIndex] = deal(NaN( length(absIdx) , 1) );
         
[R.RangeBearHead,R.SNR] = deal(NaN( length(absIdx) , 3));




% LOOP OVER RANGE CELLS

% initialize placeholder index (zero based)
pl = 0;

for rdx = 1:CS.Header.nRangeCells
    
    if ~isempty(peakIdx{rdx})
                
        % run on each range cell
        S = music_on_range_cell(CS,A,rdx,peakIdx{rdx},Vrad);
                
        % store using struct_store after some processing
        [R,pl] = store_rad_data(R,S,pl);
                
    end
end




% Compute SNR
R.SNR = get_SNR(CS,absIdx);

 
% convert to bearing to ccwE
R.RangeBearHead(:,2) = cwN2ccwE(R.RangeBearHead(:,2));


% compute lonlat:
[ekm,nkm] = rngbear2km(R.SiteOrigin,R.SiteOrigin, ...
                      R.RangeBearHead(:,1),R.RangeBearHead(:,2));
[lon,lat] = km2lonlat(R.SiteOrigin(1),R.SiteOrigin(2),ekm,nkm);

R.LonLat = [lon(:) lat(:)];


% Populate U, V
[R.U, R.V, R.Error, R.Flag] = deal(NaN(size(R.RadComp)));

% Need to compute heading ...
% * HERE *

% From loadRDLfile, unchecked really
R.U = R.RadComp .* cosd(R.RangeBearHead(:,3));
R.V = R.RadComp .* sind(R.RangeBearHead(:,3));

% Add CS file name
R.CSFileName = CS.FileName;



end

function R = get_radial_struct(flist,A)
% GET RADIAL STRUCT
%
% RADIALstruct.m plus some locally needed details


R = RADIALstruct(1);

% additional data structs
[R.TestResults, R.DOA, R.MusicParamValues, R.freqIndex, ...
                                        R.CSFileName, R.SNR] = deal([]);

% Populate expand the struct
NM = cosFileNameParts(flist{1});
R.Type = 'Meas';
R.SiteName = NM.SiteName;
R.TimeZone = 'GMT';
R.OtherMetadata.PatternUsed = A.FileName;
R.ProcessingSteps{end+1} = mfilename;
R.SiteOrigin = A.SiteOrigin;

R.README.TestResults = 'True if single brg soln';
R.README.MusicParamValues = 'Codars MUSIC parameters to determine if dual or single brg';
R.README.freqIndex = 'Frequency bin index';
R.README.DOA = '';
R.README.SNR = 'Matrix containing the SNR (col 1 is A13, 2 is A23, 3 is A33)';

% expand the struct
R(1:numel(flist)) = deal(R);


end

function  [R,pl] = store_rad_data(R,S,pl)
% STORE RAD DATA 


% append range and bearing
S.RangeBearHead = [S.Range S.Bear NaN(size(S.Bear))];

% clean up
S = rmfield(S,{'Range','Bear'});


% Insert into preallocated struct R(i)
[R,pl] = struct_store(R,pl,S,'RadComp');




end




% TESTING CODE
function test_case
% TEST CASE
%
% Compare a CSS file and the radial shorts from it with output from this
% file
% 
% Might want to build a test based on readings in Tuncer and Friedlander


% GET DATA

% Get apm
A = load_pattern_file('/m_files/test_data/music/RadialConfigs/MeasPattern.txt');

fn = '/m_files/test_data/music/CSS_KNOW_09_07_08_2220.cs4';

% Compute radials
R = radial_from_cs(fn,A);

% load comparison data
S = loadRDLFile('/m_files/test_data/music/RDLy_KNOW_2009_07_08_2220.ruv');


% MAKE PLOTS

% Show locations of solutions
figure
pws_map
hold on
h1 = plot(S.LonLat(:,1),S.LonLat(:,2),'r*');
h2 = plot(R.LonLat(:,1),R.LonLat(:,2),'bo');

legend([h1 h2],'SeaSonde','Matlab')

% Show vector map
% radial_plot(R)
% hold on
% h = radial_plot(S);
% 
% set(h,'Color','r')
% 
% pws_map


% Plot vs Bearing
figure
h1 = radial_plot_vs_bearing(S,10,1,'ro');
hold on
h2 = radial_plot_vs_bearing(R,10,1,'bo');

title('Solutions for Range Cell 10')

legend([h1 h2],'SeaSonde','Matlab')


% Plot solutions along a bearing
i = find(S.RangeBearHead(:,2) == 250);
j = find(R.RangeBearHead(:,2) == 250);

figure
h1 = plot(S.RangeBearHead(i,1),S.RadComp(i),'-r.');
hold on
h2 = plot(R.RangeBearHead(j,1),R.RadComp(j),'-b.');

title('Solutions along Bearing 250^o')

legend([h1 h2],'SeaSonde','Matlab')



% GET SHIP DATA TO TEST SINGLE/DUAL!



keyboard

% make a plot showing the peak and which range cell
if testFigs
    figure
    plot(freqs,10*log10(CS.antenna3Self(:,rdx)),'-b.'), hold on
    plot(freqs(peakIdx{rdx}),10*log10(CS.antenna3Self(peakIdx{rdx},rdx)),'g.')
    h=plot(freqs(fbin),10*log10(CS.antenna3Self(fbin,rdx)),'r*');
    title('Monopole Self Spectrum with Peaks')
    legend(h,['fbin, v=' num2str(Vrad(fbin))])
end



% USE this to check algorithm vs Example in Tony's section 4.1.
if check2
    close all
    C=[0.2162 0.0303-0.0090i 0.3170-0.0063i; ...
        0.0303+0.0090i 0.0436  -0.0091+0.0213i; ...
        0.3170+0.0063i -0.0091-0.0213i 0.5416];
    
    APM = make_ideal_pattern(225,0:5:360);
end

if check2
    disp('% ---------------------------- %')
    disp('Compare with Tonys Section 4:')
    disp('Eigen Vectors')
    eigVectors
    disp('Eigen Values')
    eigValues
end



if testFigs
    % Make figure like Tonys' figure 9
    testFig_tonysFig9(APM,DOA,singleIdx,dualIdx)
end

% Maybe run the special case for checking with Tony's paper:
if check2
    disp('% - Single Brg Case --------------%')
    disp('Magnitude^2 of the projections onto the signal and noise subspaces at 225^o are:')
    % get signal projection at 225 deg which is not normally computed:
    bdx=find(APM.BEAR==225);
    disp(['[' num2str(APM.A13R(bdx)+APM.A13I(bdx)) ' ' num2str(APM.A23R(bdx)+APM.A23I(bdx)) ' 1]'])
    A=[APM.A13R(bdx)+APM.A13I(bdx) APM.A23R(bdx)+APM.A23I(bdx) 1];
    E3=eigVectors(:,3);
    real(conj(A)*E3*conj(E3)'*A')% signal projection
    keyboard % <-- check use of ' or.'?
    1/(DOA(bdx,1)) % noise projection
    
    disp('% - Dual Brg Case --------------%')
    disp('Magnitude^2 of the projections onto the signal and noise subspaces are:')
    bdx=find(APM.BEAR==205);
    disp(['at ' num2str(APM.BEAR(bdx)) 'deg:'])
    disp(['[' num2str(APM.A13R(bdx)+APM.A13I(bdx)) ' ' num2str(APM.A23R(bdx)+APM.A23I(bdx)) ' 1]'])
    1/(DOA(bdx,2))
    bdx=find(APM.BEAR==330);
    disp(['at ' num2str(APM.BEAR(bdx)) 'deg:'])
    disp(['[' num2str(APM.A13R(bdx)+APM.A13I(bdx)) ' ' num2str(APM.A23R(bdx)+APM.A23I(bdx)) ' 1]'])
    1/(DOA(bdx,2))
end


            if ~tResult && testFigs% single brg solution case
                textBox(Pvalues,'SINGLE')
                pause, if ~check2, close all, end
            elseif tResult && testFigs % dual bearign case
                textBox(Pvalues,'DUAL')
                pause, if ~check2, close all, end
            end

            if check2
                keyboard
            end



end

function textBox(Pvalues,lab)
% add text box to figure showing music parameters
a=axis;
x=a(1)+(.1*(a(2)-a(1)));
y=a(3)+(.9*a(4)-a(3)); y2=a(3)+(.8*a(4)-a(3));
text(x,y,[ lab ' Brg Soln:'])
text(x,y2,['P=[' num2str(Pvalues,4) ']'])
return

% function Enew=complexCombine(E)
% % combine real and imag parts of each column of a matrix:
% keyboard
% Mag=sqrt( (real(E).^2) + (imag(E).^2) );
% Phas=atan2(imag(E),real(E));
% Enew=Mag.*exp(i*Phas);
%
end

function testFig_tonysFig9(APM,DOA,singleIdx,dualIdx)
% testFig_tonysFig9.M Makes a figure showing the DOA function vs bearing
% ... just like his fig 9.

figure
subplot(211)
h1=plot(APM.BEAR,10*log10(DOA(:,1)),'-b.'); hold on
plot(APM.BEAR(singleIdx),10*log10(DOA(singleIdx,1)),'g*')
title('Single Bearing DOA function')
xlabel('BEARing (deg CWN)'),ylabel('10*log10(DOA)')
text(APM.BEAR(singleIdx)+10,10*log10(DOA(singleIdx,1)), ...
    ['(' num2str(APM.BEAR(singleIdx)) ',' num2str(10*log10(DOA(singleIdx,1))) ')'])

subplot(212)
h2=plot(APM.BEAR,10*log10(DOA(:,2)),'-r.'); hold on
plot(APM.BEAR(dualIdx),10*log10(DOA(dualIdx,2)),'g*');
title('Dual Bearing DOA function')
xlabel('bearing (deg CWN)'),ylabel('10*log10(DOA)')
text(APM.BEAR(dualIdx(1))+10,10*log10(DOA(dualIdx(1),2)), ...
    ['(' num2str(APM.BEAR(dualIdx(1))) ',' num2str(10*log10(DOA(dualIdx(1),2))) ')'])
text(APM.BEAR(dualIdx(2))+10,10*log10(DOA(dualIdx(2),2)), ...
    ['(' num2str(APM.BEAR(dualIdx(2))) ',' num2str(10*log10(DOA(dualIdx(2),2))) ')'])

end

function plotMusicParams(R)
% look at properties of the music paramters (20 10 3)
%
% 1) the ratio of the largest eigenvalue to the 2nd largest is less
% than P1 (eg 20). That is (largest/2nd largest < P1).
% 2) "the ratio of the largest two signal powers to the smallest [of the
% two signal powers] is be less than P2 (usually 10)".
% 3) for the signal matrix, the ratio of the product of the diagonal
% elements to the product of the off diagonal elements must be greater
% than P3.
figure

subplot(131)
h=hist5(R.MusicParamValues(:,1),0:1:45);, set(h,'Color','b')
title('Histograms MUSIC Parameter Test Outputs')
xlabel('Eigen Value Ratios, Dual Brg if < 20)'), ylabel('Count')
hold on, a=axis; plot([20 20],[a(3:4)],'r')

subplot(132)
h=hist5(R.MusicParamValues(:,2),0:1:25);, set(h,'Color','b')
xlabel('Signal Power Ratios, Dual Brg if < 10)'), ylabel('Count')
hold on, a=axis; plot([10 10],[a(3:4)],'r')

subplot(133)
h=hist5(R.MusicParamValues(:,3),0:1:45);, set(h,'Color','b')
xlabel('Signal Matrix Ratios, Dual Brg if >3)'), ylabel('Count')
hold on, a=axis; plot([3 3],[a(3:4)],'r')
set(gcf,'units','inches','position',[0.3202    3.6691   17.5315    5.7504])

end

function plotSizer
% local version ...
set(gcf,'units','inches','Position',[2.5083    1.2141   13.3821    9.4729])
set(gca,'units','normalized','Position',[0.1300    0.1100    0.7750    0.8150])
end



% COMPARING My vs Codars Method, PLots
% move to projects or something
function checkPlotMineVsCodars(R,RADIAL)
% MAke several Plots comparing the output of my music implementation with
% that from a codar Rm file (loaded into a RADIAL structure)

% RADIAL MAPS ----
% My Music CSa:
radial_plot3('cop1', R.RangeBear(:,1), R.RangeBear(:,2), ...
    R.RadComp, [])
hold on
% COS music
radial_plot3(RADIAL.SiteName, RADIAL.RangeBearHead(:,1), ...
    RADIAL.RangeBearHead(:,2),  RADIAL.RadComp, [],'g')
axis([-121 -119 33.6 34.8])
plotSizer
ht=title(['Radials ' CSestr(RADIAL.TimeStamp) ': MyMUSIC (blk) and COS (g)']);

% PLOTS VS BEARING ----
% this Makes 3 subplots using the subfuncitons
rcells=[3 4 5];
plotVsBrg(RADIAL,R,rcells)
rcells=[23 24 25];
plotVsBrg(RADIAL,R,rcells)
rcells=[43 44 45];
plotVsBrg(RADIAL,R,rcells)
rcells=[50 51 52];
plotVsBrg(RADIAL,R,rcells)

% HISTOGRAMS OF BEARING ------
figure
subplot(211)
h=hist5((R.RangeBear(:,2)),0:360);, set(h,'Color','b')
title('Histograms of Bearing')
xlabel('Degrees CCWE'), ylabel('Count (MINE)')
subplot(212)
h=hist5(RADIAL.RangeBearHead(:,2),0:360);, set(h,'Color','r')
xlabel('Degrees CCWE'), ylabel('Count (COS''s)')
set(gcf,'units','inches','position',[2.7618    1.6411    8.1920    8.9659])

% HISTOGRAMS OF RADIAL.RadComp ------
figure
%subplot(211)
h1=hist5(R.RadComp,-90:90);, set(h1,'Color','b')
title('Histograms of Radial Component')
% xlabel('cm s^-^1'), ylabel('Count (MINE)')
%subplot(212)
h2=hist5(RADIAL.RadComp,-90:90);, set(h2,'Color','r')
xlabel('cm s^-^1'), ylabel('Count')
set(gcf,'units','inches','position',[2.7618    1.6277   14.8231    8.9792])
legend([h1 h2],'My RadComp','COS RadComp')

% Other METRIC PLOTS -------
% show vs Range Cell:
% number of common vel-range-bearing points
% number of same velcovs range cell, with diff of bearing
% number of velocities of each

% loop over range cells ...
CSaPtsCos=[];
CSaPtsMy=[];
numCommonVel=[];
numCommonVelBrg=[];

RangesCos=sort(unique(RADIAL.RangeBearHead(:,1)));
RangesMy=sort(unique(R.RangeBear(:,1)));
currentCell=unique(round([RangesCos' RangesMy']));
nRangeCells=length(currentCell);
for idx=1:nRangeCells

    % get indicies of current range cell
    i=find(round(RADIAL.RangeBearHead(:,1))==currentCell(idx));
    mdx=find(round(R.RangeBear(:,1))==currentCell(idx));

    CSaPtsCos=[CSaPtsCos length(i)];
    CSaPtsMy=[CSaPtsMy length(mdx)];

    % get number of velocity CSa points in common (rounded to 1 cm/s)
    cmnVel=intersect(round(RADIAL.RadComp(i)*2),round(R.RadComp(mdx)*2));
    % get number of velocity AND Bearign CSa points in common
    cmnVelBrg=intersect([RADIAL.RangeBearHead(i,2) round(RADIAL.RadComp(i)*2)], ...
        [R.RangeBear(mdx,2) round(R.RadComp(mdx)*2)],'rows');
    numCommonVel=[numCommonVel size(cmnVel,1)];
    numCommonVelBrg=[numCommonVelBrg  size(cmnVelBrg,1)];
end

figure
h1=plot(1:nRangeCells,CSaPtsMy,'-bo'); hold on,
h2=plot(1:nRangeCells,CSaPtsCos,'-r.');
h3=plot(1:nRangeCells,numCommonVel,'-g.');
h4=plot(1:nRangeCells,numCommonVelBrg,'-m.');
legend([h1 h2 h3 h4],'CSaPoints (my)','CSaPoints (COS)', ...
    'Common Vel (+/-1 cm s^-^1)','Com Vel&Brg')

% RANGE - VELOCITY PLOT
figure
h1=plot(R.RangeBear(:,1),R.RadComp,'bo'); hold on
h2=plot(RADIAL.RangeBearHead(:,1),RADIAL.RadComp,'r.');
legend([h1 h2],'My RadComp','COS RadComp')
xlabel('Range (km'), ylabel('cm s^-^1')

end

function plotVsBrg(RADIAL,R,rcell)
% PLOTS VS BEARING ----
% This is way cleaner ...
figure
subplot(311)
subPlotVsBrg(RADIAL,R,rcell(1))
title('Vel vs Brg: MyMUSIC (b) vs COS (r)')

subplot(312)
subPlotVsBrg(RADIAL,R,rcell(2))

subplot(313)
subPlotVsBrg(RADIAL,R,rcell(3))

xlabel('Bearing (^oCCWE)')
set(gcf,'units','inches','Position',[2.5083    1.2141   13.3821    9.4729])

end

function subPlotVsBrg(RADIAL,R,rcell)
% Run this for each subplot

% My Music CSa
rdx=find(round(R.RangeBear(:,1)*100)==round(R.rangeInc*rcell*100));
plot(R.RangeBear(rdx,2),R.RadComp(rdx),'-b.'), hold on
% add all my music radials pre- spatial mean (ie they passed parameter test
% but the were averaged to produce the final result
pdx=find(R.TestResults>0);
rdx=find(round(R.allRangeBear(pdx,1)*100)==round(R.rangeInc*rcell*100));
plot(cwN2ccwE(R.allRangeBear(pdx(rdx),2)),R.allRadComp(pdx(rdx)),'co');
myTot=length(pdx(rdx)); % p for passed test, r for this range cell

% Codar's CSa
rangeInc=min(RADIAL.RangeBearHead(:,1));
Rdx=find(round(RADIAL.RangeBearHead(:,1)*100)==round(rangeInc*rcell*100));
plot(RADIAL.RangeBearHead(Rdx,2),RADIAL.RadComp(Rdx),'-r.')
% add min/max velocity to cos CSa
Rdx=find(round(RADIAL.CS.rnge*100)==round(rangeInc*rcell*100));
plot(cwN2ccwE(RADIAL.CS.BEAR(Rdx)),RADIAL.CS.minv(Rdx),'mo')
plot(cwN2ccwE(RADIAL.CS.BEAR(Rdx)),RADIAL.CS.maxv(Rdx),'mo');
cosTot=sum(RADIAL.CS.edvc(Rdx));

axis([140 360 -90 90])
% Count up total number of veclocity solutions in this range cell, and add
% min and max velocities in range cell
text(330,70,['myTotal=' num2str(myTot)])
text(330,50,['CosTotal=' num2str(cosTot)])
text(330,30,['myRange=[' num2str(min(R.allRadComp(pdx(rdx))),3) ' ' ...
    num2str(max(R.allRadComp(pdx(rdx))),3) ']'])
text(330,10,['CosRange=[' num2str(min(RADIAL.CS.minv(Rdx)),3) ' '  ...
    num2str(max(RADIAL.CS.maxv(Rdx)),3) ']'])

ylabel(['RC ' num2str(rcell) ' cm/s'])

end


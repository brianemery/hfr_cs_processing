function APM = make_ideal_pattern(loop1Brg,bear) 
%MAKE_IDEAL_PATTERN create an ideal pattern 
% APM = make_ideal_pattern(loop1Brg,bear) 
%
% Makes an idealized version of the APM ... fields as col vectors
%
% INPUTS
% loop1Brg - the loop 1 Bearing in deg CWN ("True")
% bear     - (optional) a range of cwN bearings to cover (col vector)
%
% EXAMPLE:
% APM = make_ideal_pattern(225);

% Feb 2015 Brian M. Emery
% New version is basically a copy of Codar's, with struct output


% CODARS VERSION
%  ComputeIdealPat
% %// Compute the ideal antenna pattern
% %// The angles, Theta, are measured relative to the coastline bisector.
% %// So, 0 degrees marks the bisector and positive is toward loop 1 (clockwise)
% %// and negative is toward loop 2 (counter clockwise).
% CSSimGlobals;
% 
% %//Compute bearing angles of the pattern points
% simTh = (simFirstB:simIncrB:simLastB)';	%//'
% %//Loop 1 response
% simA1 = cos( (simTh-simLoop1B)*pi/180 ) + i*0;
% %//Loop 2 response
% simA2 = cos( (simTh-(simLoop1B-90))*pi/180 ) + i*0;
% %//Monopole response
% simA3 = ones(length(simA1),1) + i*0;
% %//Bearings relative to coast normal and in radians
% simTh = (simTh-simCoastN) * (pi/180);

% check for test case
if strcmp('--t',loop1Brg), test_case, return, end


if nargin < 2
    bear = 1:360;
end

bear = bear(:);

% Init output structure
APM = apm_struct;


% Put data into place
APM.BEAR = bear;

APM.A13R = cos( (bear-loop1Brg)*pi/180 );

APM.A23R = cos( (bear-(loop1Brg-90))*pi/180 );

% Monopole not usually used but put it in anyway
APM.A33R = ones(size(bear));

% Imaginary parts, Q factor stuff
fn = {'A13I','A23I','A33I','A13RQ','A13IQ','A23RQ','A23IQ'}; 
for i = 1:numel(fn)
    APM.(fn{i}) = zeros(size((bear))); 
end


% Mag and phase
APM = realImag2MagPhase(APM);


% loop 1 cwN:
APM.loop1Brg = loop1Brg;

% Details
APM.Type = 'Ideal Antenna Pattern';
APM.README.BEAR = 'deg cwN';
APM.Units.BEAR = 'deg cwN';

APM.ProcessingSteps = {mfilename};

end

function older_trash

% SORT OUT BEARINGS


% adjust relative to loop 1:
APM.BEAR = bear + loop1Brg;

% make anything bigger than 360 it's <360 equivalent ...
APM.BEAR(APM.BEAR >= 360) = APM.BEAR(APM.BEAR >= 360)-360;




% Subsample if necessary
if nargin == 2
   
    [c,ia,ib] = intersect(APM.BEAR,bear);
    APM = subsref_struct(APM,ia,size(APM.BEAR,1),1);
end








return

% Needs to be tested and fixed ... hack in the mean time:

APM = load_pattern_file('/m_files/test_data/IdealPattern.txt');


% adjust relative to loop 1:
APM.BEAR = APM.BEAR + loop1Brg;


% make anything bigger than 360 it's <360 equivalent ...
APM.BEAR(APM.BEAR >= 360) = APM.BEAR(APM.BEAR >= 360)-360;

% Subsample if necessary
if nargin == 2
   
    [c,ia,ib] = intersect(APM.BEAR,bear);
    APM = subsref_struct(APM,ia,size(APM.BEAR,1),1);
end

% loop 1 cwN:
APM.loop1Brg = loop1Brg;

% Details
APM.Type = 'Measured Antenna Pattern';
APM.CreateTimeStamp = now;


return


end

function put_in_test

% from an IdealPatt.txt file ...:
% blocks of idealAntennaPatt.txt are:
% loop1/3 real, realQ, imag, imagQ, 
% loop2/3 real, realQ, imag, imagQ]
% NOTES THESE ARE DEG CCW relative to loop1 bearing
brgs=[ ...
    0.   5.  10.  15.  20.  25.  30. ...
    35.  40.  45.  50.  55.  60.  65. ...
    70.  75.  80.  85.  90.  95. 100. ...
   105. 110. 115. 120. 125. 130. 135. ...
   140. 145. 150. 155. 160. 165. 170. ...
   175. 180. 185. 190. 195. 200. 205. ...
   210. 215. 220. 225. 230. 235. 240. ...
   245. 250. 255. 260. 265. 270. 275. ...
   280. 285. 290. 295. 300. 305. 310. ...
   315. 320. 325. 330. 335. 340. 345. ...
   350. 355.];
a13r=[ ...
   0.100E+01  0.996E+00  0.985E+00  0.966E+00  0.940E+00  0.906E+00  0.866E+00 ...
   0.819E+00  0.766E+00  0.707E+00  0.643E+00  0.574E+00  0.500E+00  0.423E+00 ...
   0.342E+00  0.259E+00  0.174E+00  0.872E-01 -0.271E-19 -0.872E-01 -0.174E+00 ...
  -0.259E+00 -0.342E+00 -0.423E+00 -0.500E+00 -0.574E+00 -0.643E+00 -0.707E+00 ...
  -0.766E+00 -0.819E+00 -0.866E+00 -0.906E+00 -0.940E+00 -0.966E+00 -0.985E+00 ...
  -0.996E+00 -0.100E+01 -0.996E+00 -0.985E+00 -0.966E+00 -0.940E+00 -0.906E+00 ...
  -0.866E+00 -0.819E+00 -0.766E+00 -0.707E+00 -0.643E+00 -0.574E+00 -0.500E+00 ...
  -0.423E+00 -0.342E+00 -0.259E+00 -0.174E+00 -0.872E-01  0.190E-18  0.872E-01 ...
   0.174E+00  0.259E+00  0.342E+00  0.423E+00  0.500E+00  0.574E+00  0.643E+00 ...
   0.707E+00  0.766E+00  0.819E+00  0.866E+00  0.906E+00  0.940E+00  0.966E+00 ...
   0.985E+00  0.996E+00];
% a13rQ=[
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00
a13i=[ ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00];
% a13iQ=[
%    0.389E-38  0.224E-38  0.144E-37  0.266E-28  0.168E-37  0.440E+07  0.128E-09
%    0.306E-28  0.145E+10  0.139E+11  0.159E+30  0.148E-30  0.189E+24  0.159E+30
%    0.117E-09  0.390E-38  0.404E-27  0.536E-38  0.258E+24  0.198E+14  0.606E+20
%    0.646E+24  0.239E+24  0.572E+20  0.266E-28  0.362E-38  0.159E+30  0.148E-30
%    0.159E+24  0.175E-23  0.440E+07 -0.171E+00  0.107E+15  0.160E+07  0.445E+20
%    0.266E-28  0.435E-38  0.159E+30  0.148E-30  0.159E+24  0.265E-20  0.241E+07
%   -0.189E+39  0.159E+30  0.145E+10  0.139E+11  0.644E-23  0.932E+09 -0.469E-31
%    0.185E+26  0.756E+29  0.281E-40   NAN(072)  0.279E-37  0.749E-39  0.127E-09
%    0.442E-34  0.110E-10  0.202E-18  0.159E+30  0.202E-18  0.217E-18  0.548E+20
%    0.128E-09 -0.263E+38  0.407E-24  0.443E+07  0.222E+24  0.113E-38  0.187E+24
%   -0.826E+10  0.111E-38
a23r=[ ...
   0.000E+00  0.872E-01  0.174E+00  0.259E+00  0.342E+00  0.423E+00  0.500E+00 ...
   0.574E+00  0.643E+00  0.707E+00  0.766E+00  0.819E+00  0.866E+00  0.906E+00 ...
   0.940E+00  0.966E+00  0.985E+00  0.996E+00  0.100E+01  0.996E+00  0.985E+00 ...
   0.966E+00  0.940E+00  0.906E+00  0.866E+00  0.819E+00  0.766E+00  0.707E+00 ...
   0.643E+00  0.574E+00  0.500E+00  0.423E+00  0.342E+00  0.259E+00  0.174E+00 ...
   0.872E-01 -0.542E-19 -0.872E-01 -0.174E+00 -0.259E+00 -0.342E+00 -0.423E+00 ...
  -0.500E+00 -0.574E+00 -0.643E+00 -0.707E+00 -0.766E+00 -0.819E+00 -0.866E+00 ...
  -0.906E+00 -0.940E+00 -0.966E+00 -0.985E+00 -0.996E+00 -0.100E+01 -0.996E+00 ...
  -0.985E+00 -0.966E+00 -0.940E+00 -0.906E+00 -0.866E+00 -0.819E+00 -0.766E+00 ...
  -0.707E+00 -0.643E+00 -0.574E+00 -0.500E+00 -0.423E+00 -0.342E+00 -0.259E+00 ...
  -0.174E+00 -0.872E-01];
% a23rQ=[
%    0.159E+30  0.962E-40  0.123E-38  0.482E+13  0.159E+30  0.148E-30  0.196E+24
%    0.219E+12  0.159E+30  0.421E+07  0.159E+30  0.148E-30  0.554E+24  0.266E-28
%    0.123E-38  0.219E+12  0.119E-18  0.786E-20  0.119E-18  0.509E-28  0.932E+09
%   -0.117E-31  0.241E+30  0.182E+29  0.281E-40   NAN(072)  0.322E-19   NAN(112)
%    0.110E-10  0.110E-10  0.110E-10  0.110E-10  0.110E-10  0.104E+35  0.202E-18
%    0.176E-18  0.119E+06  0.202E-18  0.132E-10   NAN(114)  0.110E-10  0.195E+07
%    0.127E-09 -0.218E+39  0.781E-20  0.266E-28  0.371E-38  0.260E+24  0.250E-30
%    0.747E-39  0.113E-38  0.250E-30  0.747E-39  0.150E-38 -0.500E+06 -0.504E+06
%   -0.168E+29  0.988E-37  0.306E-28  0.145E+10  0.159E-09  0.159E+30  0.145E+10
%    0.241E+07  0.224E-38  0.165E-35  0.157E+10   NAN(000)  0.579E-38   NAN(096)
%    0.726E-42   NAN(012)
a23i=[ ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00 ...
   0.000E+00  0.000E+00];
% a23iQ=[
%    0.140E-44  0.918E-40  0.000E+00  0.000E+00  0.000E+00  0.552E-41  0.552E-41
%    0.918E-40  0.640E-39  0.214E-37  0.134E-39  0.344E-41  0.420E-44  0.420E-44
%    0.151E-38  0.214E-37  0.409E+01  0.408E+01  0.408E+01  0.408E+01  0.000E+00
%     NAN(255)   NAN(000)  0.420E-44  0.311E-39  0.311E-39  0.436E+01  0.317E-39
%    0.312E-39  0.140E-44  0.918E-40  0.000E+00  0.000E+00   NAN(000)  0.300E-41
%    0.300E-41  0.918E-40  0.640E-39  0.214E-37  0.169E-38  0.344E-41  0.215E-37
%    0.000E+00  0.000E+00  0.276E-39  0.165E-38  0.276E-39  0.165E-38  0.276E-39
%    0.165E-38  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.918E-40  0.000E+00  0.000E+00  0.000E+00  0.000E+00
%    0.000E+00  0.000E+00  0.000E+00  0.841E-44  0.168E-43  0.266E-43  0.350E-43
%    0.448E-43  0.392E-43

% Convert brgs, which are CCW relative to loop1, to cwN 
brgsCWN = loop1Brg-brgs;

% make anything bigger than 360 it's <360 equivalent ...
brgsCWN(brgsCWN>=360)=brgsCWN(brgsCWN>=360)-360;

% ... and anything less than 0 likewise:
brgsCWN(brgsCWN<0)=brgsCWN(brgsCWN<0)+360;


% now make it so that it only covers the same range as the input apm, and
% scale same size as input magnitudes
if nargin == 2
[c,n,ib]=intersect(brgsCWN,bear); clear c ib
else
    n = 1:length(brgsCWN);
end
% a13m=a13m(i);%.*max(A1);
% a23m=a23m(i);%.*max(A2);
% phase1=phase1(i);
% phase2=phase2(i);

% put it all into the output structure
APM = apm_struct;
APM.BEAR = brgsCWN(n);
APM.A13R = a13r(n);
APM.A13I = a13i(n);
APM.A23R = a23r(n);
APM.A23I = a23i(n);
APM.A13M = sqrt( (a13r(n).^2) + (a13i(n).^2) );
APM.A13P = atan2(a13i(n),a13r(n)).*180/pi;
APM.A23M = sqrt( (a23r(n).^2) + (a23i(n).^2) );
APM.A23P = atan2(a23i(n),a23r(n)).*180/pi;
APM.loop1Brg =loop1Brg;
APM.README.BEAR = 'deg cwN';

% Meta data
APM.Type = 'Ideal Antenna Pattern';
APM.CreateTimeStamp = now;
APM.README.CreatedWith = mfilename;

keyboard

return

end

function older_test

%keyboard

% Code Checking Code:
APM=makeIdealPatt(239,[0:5:360]-1);



makeAPMplot(APM)
title('From makeIdealPattern.m')
%out=[APM.bear' APM.a13r' APM.a13i' APM.a23r' APM.a23i' APM.a13m' APM.a13p'
%APM.a23m' APM.a23p'] 

% make comparison figure 
file='/Data/SSD APM/PATT_ssd1_ideal_239deg.xlp/SEAS_ssd1_2008_02_26_0704.patt';
APM=loadPatternFile(file);
makeAPMplot(APM)
title('From SEAS_ideal.patt file')

% *Checked vs. Tony's paper and an idealPatt.txt file, as well as an ideal
% SEAS.patt file. COS may have an error, in that their output is off by one
% degree in the SEAS files
% 10jun08 -BE

% example:
APM=makeIdealPatt(225,0:5:360);
bdx=find(APM.bear==330);
disp(['[' num2str(APM.a13r(bdx)+APM.a13i(bdx)) ' ' num2str(APM.a23r(bdx)+APM.a23i(bdx)) ' 1]'])
bdx=find(APM.bear==205);
disp(['[' num2str(APM.a13r(bdx)+APM.a13i(bdx)) ' ' num2str(APM.a23r(bdx)+APM.a23i(bdx)) ' 1]'])
bdx=find(APM.bear==225);
disp(['[' num2str(APM.a13r(bdx)+APM.a13i(bdx)) ' ' num2str(APM.a23r(bdx)+APM.a23i(bdx)) ' 1]'])
end

function test_case
% TEST CASE


% COMPARE WITH IDEAL
% this uses both SEAS and IdealPattern formats, with Loop 1 at 249 deg
wd = '/m_files/test_data/make_ideal_pattern/';

I = load_pattern_file([wd 'IdealPattern_XXXX_2015_02_26_1648.txt']);
S = load_pattern_file([wd 'SEAS_XXXX_2015_02_26_1648.txt']);


keyboard

A = make_ideal_pattern(0,1:360);




% -- PREVIOUS --

% COMAPRE WITH IDEAL
% Comapare it with an ideal pattern file
I = load_pattern_file('/m_files/test_data/make_ideal_pattern/IdealPattern.txt');
H = plot_apm_polar(I);

% Load, and over plot
A = make_ideal_pattern(0,1:360);
LS = line_style_groups;
[LS(1:2).Marker] = deal('*');

H = plot_apm_polar(A,LS,H);


% Sort the file data
[y,i] = sort(I.BEAR);
II = subsref_struct(I,i,length(I.BEAR),1);

fn = {'A13R','A13I','A23R','A23I','A13M','A13P','A23M','A23P','A33R','A33I'};

for i = 1:numel(fn)
   
    rslt = isequal(round(A.(fn{i})*1e3)/1e3,round(II.(fn{i})*1e3)/1e3);
    
    disp([' test ' fn{i} ' ... isequal = ' num2str(rslt) ])
   
    keyboard
end

keyboard

% COMPARE WITH MEAS

% Compare it with a measured pattern file
M = load_pattern_file('/m_files/test_data/make_ideal_pattern/MeasPattern.txt');


figure
H = plot_apm_polar(M);


% Load, and over plot
A = make_ideal_pattern(M.loop1Brg,M.BEAR);
LS = line_style_groups;
LS(1).Marker = '*';
LS(2).Marker = '^';

H = plot_apm_polar(A,LS,H);



end

function APM = load_pattern_file(file)
%LOAD_PATTERN_FILE load several kinds of APM file (txt, SEAS, loop)
% APM = load_pattern_file(file) 
%
% Loads SEAS and MeasPattern.txt type antenna patter
% files into a structure.
%
% see also: evaluate_apm.m, and File_AntennaPattern.pdf
% and ctfReader.m

% Copyright (C) 2009-2010 Brian M. Emery
% June 2009

% TO DO
% get some units into the mix, add other file support
%
% check that fields are columns for all types
% 
% I dont think this works with IdealPattern files?
%
% SORT THE APM IN BEAERING
% [th,ix] = sort(APM.BEAR);
% w
% APM
% APM = subsref_struct(APM,ix,size(APM.BEAR,1),1)

% check for test case
if strcmp('--t',file), test_case, return, end


% Recurse for n cell inputs
if iscell(file)
    for i = 1:numel(file)
        APM(i) = load_pattern_file(file{i});
    end
    return
end

% get file name parts
[pathstr,fname,ext] = fileparts(file); 


% load the compatible file types
if strncmp('SEAS',fname,4) 
    APM = load_seas_file(file);
    
elseif any(strncmp({'Meas','Idea'},fname,4))  && strcmp('.txt',ext)
    APM = load_meas_patt_file(file);
    
elseif strncmp('LOOP',fname,4)
    
    % Get the loop file ... crude and needs work
    S = ctfReader(file);
    APM = loop_to_apm(S);
    
else
    disp([fname ' is not a supported type']), APM=[]; return
end

eval(['APM.FileName=''' file ''';'])

end


function APM = load_seas_file(file)
% LOAD SEAS FILE.M
% out=load_seas_file(file)
% Load SEAS pattern file data in and get it organized (SEAS_*.patt FILES CASE)
% SEAS_*.patt files contain:
% Bearing    A13        A13        A23        A23     Flag     A13      A13      A23      A23 
%  DegCW     real       imag       real       imag             mag      deg      mag      deg 
%
% NOTE That the bearings in SEAS files seem to be CWN relative to the Loop
% 1 direction ....

% Modifications 15 Oct 2010: updated methods to use ctfReader.m

% Use CTFREADER, produces nearly same structure as apm_struct.m
APM = ctfReader(file);

% Get loop 1 bearing. The bearing data in BEAR is relative to this number,
% where the data in BNCW is relative to north.
[idx,nm,val]=getNameValuePair('AntennaBearing',APM.HdrNames,APM.HdrValues);
strs=strparser(val{:});
APM.loop1Brg = str2num(strs(1,:)); 

% Standardize BEAR data
if isfield(APM,'BNCW')
    APM.BEAR = APM.BNCW;
    APM = rmfield(APM,'BNCW');
      
else
    % older SEAS files dont have BNCW field
    APM.BEAR = APM.BEAR + APM.loop1Brg;

    % make anything bigger than 360 it's <360 equivalent
    APM.BEAR(APM.BEAR >= 360)=APM.BEAR(APM.BEAR >= 360)-360;

end

% Meta
APM.README.BEAR_Units = 'degCWN';
APM.README.loop1Brg_Units = 'degCWN';
APM.ProcessingSteps{end+1} = mfilename;





end

function [loop1Brg,loop1BrgUnits] = getLoop1Brg(file)
% GET LOOP 1 BEARING
% Get the loop 1 direction in deg CWN.
% This could probably be generalized to work with a lot of COS data.

[names, values] = getRDLHeader(file);
[idx,nm,val]=getNameValuePair('AntennaBearing',names,values);
strs=strparser(val{:});
loop1Brg = str2num(strs(1,:)); 
loop1BrgUnits = val;

end

function APM = load_meas_patt_file(file)
% APM=load_meas_patt_file(file)
% Loads MeasuredPatt.txt file into an APM structure
%
% measuredPattern.txt file:
% Bearings CCW from Header loop1 bearing (CWN)
% Contents:
% line1: number of bearings
% Bearing deg CCW from loop 1
% Loop1to3 real
% Loop1to3 real std dev of measurments (quality factor)
% Loop1to3 imaginary
% Loop1to3 imaginary quality factor
% Loop2to3 real
% Loop2to3 real std dev of measurments (quality factor)
% Loop2to3 imaginary
% Loop2to3 imaginary quality factor
%
% Bearings in the files are CCW relative to loop1 bearing, which SeaSonde
% software gets from the Header.txt file.
%
% Based on suggestion by Simone Cosoli 2008, using COS's documentation, as
% implemented by Brian Emery

% Copyright (C) 2011  Brian Emery

% blocks of .txt are DEG CCW relative to loop1 bearing

% INIT APM struct
APM = apm_struct;

% Load file in (based on Simone's idea)
fid = fopen(file,'r');
nBearing = str2num(fgetl(fid));

bearAng = fscanf(fid,'%f',[nBearing,1]);
APM.A13R   = fscanf(fid,'%f',[nBearing,1]);
APM.A13RQ  = fscanf(fid,'%f',[nBearing,1]);
APM.A13I   = fscanf(fid,'%f',[nBearing,1]);
APM.A13IQ  = fscanf(fid,'%f',[nBearing,1]);
APM.A23R   = fscanf(fid,'%f',[nBearing,1]);
APM.A23RQ  = fscanf(fid,'%f',[nBearing,1]);
APM.A23I   = fscanf(fid,'%f',[nBearing,1]);
APM.A23IQ  = fscanf(fid,'%f',[nBearing,1]);


% Insert these ... sometimes I need them
APM.A33R = ones(size(APM.A13R));
APM.A33I = zeros(size(APM.A13R));


% Get the footer contents and parse out stuff 
footer={};
 while 1
     tline = fgetl(fid);
     
     if ~ischar(tline)
         break
     end
     
     if ~isempty(tline)
         footer{end+1} = tline;
     end
 end
fclose(fid);


% Get variable 'names' and 'values' from footer
% (David's code (from loadRDLfile))
[names,values] = deal(cell(length(footer),1));
for k = 1:length(footer)
    % This way is considerably more efficient than using strtok
    ii = find(footer{k} == '!');
    values{k,1} = strtrim(footer{k}(1:ii-1)); 
    names{k,1}  = strtrim(footer{k}(ii+1:end)); % Removie initial '!'
end

% Get the loop 1 bearing
idx=strmatch('Antenna Bearing',names);
APM.loop1Brg = str2double(values{idx});

% NOTE bearAng is DEG CCW relative to loop1 bearing
APM.BEAR = ccwE2cwN(bearAng+cwN2ccwE(APM.loop1Brg));

% Compute Mag, Phase
APM = realImag2MagPhase(APM);

% % Get Magnitudes
% APM.A13M = sqrt( (APM.A13R.^2) + (APM.A13I.^2) );
% APM.A23M = sqrt( (APM.A23R.^2) + (APM.A23I.^2) );
% 
% % Get Phases
% APM.A13P = atan2(APM.A13I,APM.A13R).*180/pi;
% APM.A23P = atan2(APM.A23I,APM.A23R).*180/pi;

% Add to README:
APM.README.BEAR_Units = 'degCWN';
APM.README.loop1Brg_Units = 'degCWN';
APM.Units.BEAR = 'degCWN';

% Add footer 
APM.Footer.Names  = names;
APM.Footer.Values = values;
APM.Footer.Orig = footer;

% Add file name
APM.FileName = file;


% Take care of some details
APM.ProcessingSteps{end+1} = mfilename;
APM.ProcessingSteps{end+1} = 'realImag2MagPhase';

idx = strmatch('Site Code',APM.Footer.Names);
if ~isempty(idx)
    APM.SiteName = APM.Footer.Values{idx};
end

idx = strmatch('Site Lat Lon',APM.Footer.Names);
if ~isempty(idx)
    APM.SiteOrigin = str2num(APM.Footer.Values{idx});
    APM.SiteOrigin = [APM.SiteOrigin(2) APM.SiteOrigin(1)]; 
end

idx = strmatch('Date Year Mo Day Hr Mn Sec',APM.Footer.Names);
if ~isempty(idx)
    APM.CreateTimeStamp = datestr(datenum(APM.Footer.Values(idx),'yyyy mm dd HH MM SS')) ;
end

%disp(['Read ' file])


return
 

end




function test_case
% TEST CASE
%
% Basic functionality test ...

% TO DO
% test the bearing adjustment on older seas files, 

% COMPARE WITH IDEAL
% this uses both SEAS and IdealPattern formats, with Loop 1 at 249 deg
wd = '/m_files/test_data/make_ideal_pattern/';

I = load_pattern_file([wd 'IdealPattern_XXXX_2015_02_26_1648.txt']);
S = load_pattern_file([wd 'SEAS_XXXX_2015_02_26_1648.patt']);

H = plot_apm(I);

figure

H = plot_apm(S);



keyboard

% investigate a very small diff
for i = 1:10, figure(1), pause, figure(2), pause, end


end

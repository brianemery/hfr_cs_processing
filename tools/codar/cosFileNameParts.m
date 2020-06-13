function  NM = cosFileNameParts(fileName)
% COS FILE NAME PARTS - get info from Codar file names
% Takes a compatible file name (cell or char) and parses it into parts, 
% outputing these into a structure with fields, for example:
% namePartStruct = cosFileNameParts(fileName)
% 
% namePartStruct = 
% 
%             Type: 'RDLm'
%         SiteName: 'CLHB'
%        TimeStamp: 7.3334e+05
%          yrmoStr: '0710'
%     wholeDateStr: '2007_10_29_210013'
%       oldDateStr: '07_10_29_2100'
%
% Works on one file at a time (char inputs) or multi-element cell arrays.
%
% See parseRDLFileName.m which seems to do similar things - might be able
% to combine these ...
% 
% or: nm = dir([csqDataDir '/CSQ_*']);
% csqNames=char({nm.name}); 
% csqTimes=datenum(csqNames(:,end-17:end-3),'yy_mm_dd_HHMMSS');
%
% Run a test case:
% NM = cosFileNameParts('--t')

% TO DO
% June 2015: generalize like this fnames_to_times.m? (datenum slow?)
%
% Could be a lot faster, see:
% csqListTimes = fnames_to_times(csqListNames,['CSQ_' AIS.SiteName '_'],'yy_mm_dd_HHMMSS');
%
% handle files like LOOP_test.txt gracefully
% allow cell inputs, vectorize fileparts_my.m
% maybe convert everything to number and matricies that way?

% Copyright (C) 2008-2010 Brian M. Emery
% 14 April 2008 
% 02 June 2010 added test case
% 08 Augu 2010 added support for cell inputs
% 02 Dec 2010 Added use of datestrq.m, improved for speed
%             Added seconds to wholeDateStr


% Optionally run the test case
if strcmp(fileName,'--t')
    NM = test_case; return
end

% --------------------------------------------------------- 
%  RECURSE FOR MULTI-ELEM INPUTS
%---------------------------------------------------------- 

% Iteratively accept multi0-element cell inputs. Probably not 
% the fastest, but it's quick to code ...
if iscell(fileName) 
    
    % Initialize the output
    NM = init_struct(numel(fileName));
    
    for i = 1:numel(fileName)
        NM(i) = cosFileNameParts(fileName{i});
    end 
    return
end


% --------------------------------------------------------- 
%  INITIALIZE 
%---------------------------------------------------------- 

% initialize structure (use subfunction)
NM = init_struct;

% identify the (potentialy) supported file name types, 
% 3 and 4 character long:
types={'STAT','RDLi','RDLm','Radz','Rads','DIAG','LOOP',...
       'TRAK','P_Rad','RDLy','ELTm','ELTi','TOTL','SEAS'}; 
types3={'CSA','Rng','Lvl','CSS','CSQ'};
    

% --------------------------------------------------------- 
%  PARSE FILE NAME 
%---------------------------------------------------------- 

% if input is cell, convert to char
if iscell(fileName)
    fileName=char(fileName);
end

% get important parts
[pathstr,fileName,ext] = fileparts(fileName);

% figure out the current file type, with error catching
NM.Type=[types{find(strncmp(fileName,types,4))}, types3{find(strncmp(fileName,types3,3))}];
if isempty(NM.Type)
    %disp('cosFileNameParts inputs ='), fileName
    disp('fileType not supported by cosFileNameParts')
    NM = []; return
end

% --------------------------------------------------------- 
%  PUT INTO A STRUCTURE 
%---------------------------------------------------------- 
% Convert name into time and site name
switch NM.Type

    
    case {'RDLi','RDLm','RDLy','ELTm','ELTi','TRAK','TOTL','SEAS'}
        stime=datenum(fileName(end-16:end),'yyyy_mm_dd_HHMM'); 
        NM.SiteName=fileName(6:9);
               
    case {'Rng','Lvl'}
        stime=datenum(fileName(10:26),'yyyy_mm_dd_HHMMSS');
        NM.SiteName=fileName(5:8);
        
    case {'CSA','CSS','Radz','Rads'} 
        stime=datenum(fileName(10:22),'yy_mm_dd_HHMM');
        NM.SiteName=fileName(5:8);

    case {'CSQ'} 
        stime=datenum(fileName(10:24),'yy_mm_dd_HHMMSS');
        NM.SiteName=fileName(5:8);

    case {'P_Rad'} 
        stime=datenum(fileName(12:24),'yy_mm_dd_HHMM');
        NM.SiteName=fileName(7:10);

    case {'STAT'}
        stime=datenum(fileName(11:20),'yyyy_mm_dd');
        NM.SiteName=fileName(6:9);
        
    case {'DIAG'}
        stime=datenum(fileName(11:18),'yyyymmdd');
        NM.SiteName=fileName(6:9);

    case {'LOOP'}
        try stime=datenum(fileName(11:23),'yymmdd_HHMMSS');
        catch
            stime = datenum(fileName(11:16),'yymmdd');
        end
        NM.SiteName=fileName(6:9);

    otherwise
        disp([fileName ' is not a supported fileType'])
    keyboard
end

% create outputs: 
NM.TimeStamp=stime;

% datestr is slow, so do it once using optimized code:
NM.wholeDateStr = datestrq(stime);%,'yyyy_mm_dd_HHMM'); % eg 2007_10_29_2100

% then use the parts
NM.yrmoStr = NM.wholeDateStr([3 4 6 7]); %'yymm'); % eg '0710'
NM.oldDateStr = NM.wholeDateStr(3:end-2); %'yy_mm_dd_HHMM'); % eg 07_10_29_2100
%NM.wholeDateStr = NM.wholeDateStr(1:end-2); % eg 2007_10_29_2100


% % datestr version
% NM.wholeDateStr=datestr(stime,'yyyy_mm_dd_HHMM'); % eg 2007_10_29_2100
% NM.yrmoStr = datestr(stime,'yymm'); % eg '0710'
% NM.wholeDateStr=datestr(stime,'yyyy_mm_dd_HHMM'); % eg 2007_10_29_2100
% NM.oldDateStr=datestr(stime,'yy_mm_dd_HHMMSS'); % eg 07_10_29_2100

end
% --------------------------------------------------------- 
function NM = test_case
% TEST CODE 

fn={'STAT_Rfg1_2007_07_08.rdt', ...
    'LOOP_cop1_100723_211027.loop', ...
    'RDLm_Rfg1_2008_08_29_1900.ruv', ...
    'TOTL_PWSS_2009_07_15_2200.tuv', ...
    'Rng_Rfg1_2007_05_14_042228.rs', ...
    'CSS_Rfg1_07_03_28_1420.cs4', ...
    'RadzRfg1_07_08_05_0800.rv', ...
    'P_Rad_sci1_08_05_03_0100.mat', ...
    'DIAG_cop1_20091031', ...
    'SEAS_ssd1_2008_02_26_0704'};

for i=1:numel(fn), fn(i), NM = cosFileNameParts(fn(i)), pause, end

% TEST CELL INPUTS
disp('testing cell inputs')
NM = cosFileNameParts(fn);

keyboard

% Second test for optimizing for speed (109.4 seconds before, 47.0 of it on datestr):
load /Volumes/CodarData/Data/testData/cosFileNameParts/csqNames.mat csqNames
csqNames = csqNames(1:3075);
profile clear
profile on
NM = cosFileNameParts(csqNames);
profile report % latest was ~21 seconds

end
% ---------------------------------------------------------
function NM = init_struct(i)


NM=struct('Type',[],'SiteName',[],'TimeStamp',[], ...
    'yrmoStr',[],'wholeDateStr',[],'oldDateStr',[]);

if nargin ==1
    NM(1:i) = NM;
end

end

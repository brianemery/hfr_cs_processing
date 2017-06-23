function [uhf,vhf,thf,grdOut,Ct,angd,ehf] = get_total_data(t_min,t_max,lon,lat,drive)  
% GET TOTAL DATA 
% [U,V,stime,gridd,Count,AngD,Err] = get_total_data(t_min,t_max,in_lon,in_lat);
%
% Extracts CODAR total velocity data at specified times and
% lon lats. Uses unix code to download data from euler.msi.ucsb.edu if it's
% not found locally (only works on a mac or unix machine).  By default,
% looks for data in /data/totals/ or c:\data\totals\, but will accept a 5th
% input with a different folder path.
%
% BACKWARD COMPATIBILITY NOTES
% Regarding the AngDiff variable, the new total data (after Jan 2011),
% processed with the HFRprogs toolbox, has several estimates of GDOP that
% should be more robust than the AngDiff. These error estimates are based
% on some pretty solid theory. The screening is applied to the TUV struct
% using a function called cleanTotals in the subfunction regular_load.
% Here, the new variables are translated into the old, and screening is
% applied. The AngDiff variable is created and set all to 90 for backward
% compatibility. 
% 
% 
% INPUTS
% lon lat inputs are optional (gets whole data set in this case)
% 
% EXAMPLE
%
% % SB Channel area
% in_lon = [-121.0 -119.0]; % longitude boundary
% in_lat = [33.8 35.0];   % latitude boundary
% st = datenum(2004,8,16,17,0,0):1/24:datenum(2004,8,19,16,0,0);
% total_dir = '/data/totals/';
% 
% [U,V,stime,gridd,Count,AngD,Err] = ...
%           get_total_data(st(1),st(end),in_lon,in_lat,total_dir);
%
% % Whole grid, on euler
% total_dir = '/Users/codar/Sites/tmp/totals/';
% in_lon = [-122.0 -118.0]; % longitude boundary
% in_lat = [33 35.5];   % latitude boundary
% st = datenum(2011,11,16,17,0,0):1/24:datenum(2011,11,19,16,0,0);
% 
% [U,V,stime,gridd,Count,AngD,Err] = ...
%           get_total_data(st(1),st(end),in_lon,in_lat,total_dir);
 

% Copyright (C) 2013 Brian  Emery

% TO DO
% better check sum - add up ~isnans ??
% HFRP version
% Structureize the loads
% use subsref_struct
% test


% SORT OUT INPUTS

% Optional test case
if strcmp(t_min,'--t')
    test_case; [uhf,vhf,thf,grdOut,Ct,angd,ehf] = deal([]); return
end

% Define the output grid based on times
grdOut = define_grid([t_min t_max]);

% Convert input lon lat to grid points
if nargin > 2
    
    ix = find(grdOut(:,1) > lon(1) & grdOut(:,1) < lon(2) & ...
              grdOut(:,2) > lat(1) & grdOut(:,2) < lat(2));
              
    % modify based on resulting grid index
    if ~isempty(ix), grdOut = grdOut(ix,:); end
    
end

% create the output time array
st = (round(t_min*24)/24):(1/24):(round(t_max*24)/24);

% compute total hours of output data
hrs = length(st);

% create the output arrays (NaN's)
[uhf,vhf,Ct,angd,ehf]= deal(NaN(size(grdOut,1),hrs));
thf=NaN*ones(1,hrs);



% GET DATA
% NO references to grdndx below, using intersect.m to get like pts 
% between grdOut and loaded gridd

% get full path list of monthly files
if nargin < 5
    totNames = times2totnames([t_min t_max]); 
else
    totNames = times2totnames([t_min t_max],drive);
end


% load the total data and save the time points you want
checkNum=0; 



% LOAD DATA

for j = 1:numel(totNames);
    
    
    % look for it in default directory
    
    if exist(totNames{j},'file')        
        
        % load it
        S = regular_load(totNames{j});
 
       
    elseif isunix && ~exist(totNames{j},'file') 
        
        % load or download and load
        S = unix_load(totNames{j});      
 
         
    elseif ~exist(totNames{j},'file') 
        
        % likely on a pc and file not found (clearer instructions needed)
        disp('file not found ...')
        disp('Grab it manually from http://euler.msi.ucsb.edu/~codar/tmp/totals/')
        return
        
    end
    
    
        

    % Get matching times
    [~,iOut,iTot] = intersect_hour_round(st,S.serialtime);
    
    checkNum = checkNum+length(iTot); 
        
    % get the matching grid points (grid pts must be same to thousanth decimal place)
    [~,iGrd,iGrdOut] = intersect(round(1e4*S.gridd),round(1e4*grdOut),'rows');
    

    % Detect and solve problems with grid intersect. Due to changes in
    % lonlat2km algorithms, and possibly due to single/double conversions,
    % the grids can be as much as 2 meters apart, causeing intersect.m to
    % not work as desired. In this case, convert to km, and round to
    % nearest 10 m
    if length(grdOut) ~= length(iGrdOut)       
        [overLap,iGrd,iGrdOut] = use_km_grid(S.gridd,grdOut);
    end
    
    
    angd(iGrdOut,iOut) = S.AngDiff(iGrd,iTot);
     uhf(iGrdOut,iOut) = S.U(iGrd,iTot);
     vhf(iGrdOut,iOut) = S.V(iGrd,iTot);
      Ct(iGrdOut,iOut) = S.Count(iGrd,iTot);
     ehf(iGrdOut,iOut) = S.Err(iGrd,iTot);
             thf(iOut) = S.serialtime(iTot);
    
    
end

% output some info to make sure the intersect line is
% working properly (this assumes st starts after sept97)
if length(st)>checkNum
    disp('WARNING!! possibly missing data')
    disp('Amount of matching data (in time) is less than requested')
    disp('Check the timeseries')
    keyboard
elseif length(st)<checkNum
    disp('WARNING!! intersect.m might not be working in get_total_data.m?')
    disp('More matching data than the number of requested times ')
    disp('check that input timeseries contains unique elements')
    keyboard
end




end

function grdOut = define_grid(st)
% DEFINE GRID
% define the grid used by the total data as the 
% array has changed over time.


% just define the max size grid and let the user decide if they want a
% smaller one?
grdOut = sbgrid(6);



% % this is the 'define specGrd' case
% if st(end) < datenum(2003,04,01)
%     
%     opt = 2;  % use 5 site grid
%     
% elseif st(1) >= datenum(2003,04,01)
%     
%     opt = 4; % 3;  % use orig. grid with east channel
% 
% % elseif st(end) >= datenum(2003,04,01) & st(1) < datenum(2003,04,01)
% % 
% %     opt = 4;  %large composite grid, appending east
% % 
% % elseif st(end) > datenum(2008,5,1)
% % 
% %     opt = 4;  %large composite grid, appending east
% 
% end
% 
%     
% 
% % % define grid 
% % if opt < 4
% 
%     grdOut = sbgrid(opt);
% 
% % elseif opt == 4
% % 
% %     % this should be done in sbgrid.m
% %     grdOut1 = sbgrid(2);
% %     grdOut2 = sbgrid(3);
% %     grdOut = [grdOut1; grdOut2(1552:end,:);]; 
% % 
% % end


end

function [overLap,iGrd,iGrdOut] = use_km_grid(gridd,grdOut)
% USE KM GRID - convert to km to match up grid points

% Get sbgrid central location
loc = codar_sites({'central'});

% conver to km
[ekm,  nkm ] = lonlat2km(loc.central(1),loc.central(2), gridd(:,1), gridd(:,2));
[ekm2, nkm2] = lonlat2km(loc.central(1),loc.central(2),grdOut(:,1),grdOut(:,2));

% now do matching to nearest 10 m
[overLap,iGrd,iGrdOut]=intersect(round(10*[ekm nkm]),round(10*[ekm2 nkm2]),'rows');

disp('used km grid')
end

function S = unix_load(totstr)
% UNIX LOAD

% get the file name only
[~,fname,~] = fileparts(totstr);
    
% look for it in downloads
if exist(['~/Downloads/' fname '.mat'],'file')
    
    S = regular_load(['~/Downloads/' fname '.mat']);
        
% otherwise try to download it
else
    
    % inform the user
    disp([totstr ' not found ...'])
    disp(['downloading from web to ~/Downloads/'])
    
    
    % get pwd
    cdir = pwd; 

    % cd to the downloads directory 
     cd('~/Downloads/')

    % download it            
    [r,s] = system(['curl -O http://euler.msi.ucsb.edu/~codar/tmp/totals/' fname '.mat'],'-echo');
    
    % go back to orig directory
    cd(cdir)
   
    % now load it
    S = regular_load(['~/Downloads/' fname '.mat']);
    
      
    
end
        
end

function A = regular_load(totstr)
% REGULAR LOAD
% 
% load command with messaging and details
%
%
% From:
% Kaplan, David M., John Largier, and Louis W. Botsford.
% "HF radar observations of surface circulation off Bodega Bay
% (northern California, USA)." Journal of Geophysical Research:
% Oceans (1978?2012) 110.C10 (2005).
%
% "Total currents were rejected if they were greater than 1.0 m/s or they
% had a geometric dilution of precision (GDOP) [Wells and Beck, 1987;
% Chapman and Graber, 1997] greater than 2. We also used a second
% estimation of the error in a total current that is similar to the mapping
% error but takes into account the measured variability in the radial
% currents [Lipa, 2003]. If this error vector was greater than 0.18 m/s in
% magnitude then the current estimate for that grid point and time period
% was omitted. Given that our CODAR radial currents have a noise level of
% order 0.06? 0.10 m/s


disp([totstr ' ... loaded']);

S = load(totstr);

if isfield(S,'TUV')
    
    
    % unpack TUV struct
    A = S.TUV;

    % created needed fields (AngDiff, Count, Err, serialtime, gridd, names)
    A.AngDiff = 90.*ones(size(A.U));
    A.Count = A.OtherMatrixVars.makeTotals_TotalsNumRads;
    A.Err = A.ErrorEstimates(1).TotalErrors; % GDOP 
    A.serialtime = A.TimeStamp;
    A.gridd = A.LonLat;
    A.names = {};
    
    % Apply error screening statistics
    % From:
    % Kaplan, David M., John Largier, and Louis W. Botsford. 
    % "HF radar observations of surface circulation off Bodega Bay 
    % (northern California, USA)." Journal of Geophysical Research: 
    % Oceans (1978?2012) 110.C10 (2005).
    A = cleanTotals(A,100,{ 'GDOP','TotalErrors',2 },{ 'GDOPRadErr','TotalErrors',18 });
    
    % clear redundant TUV fields?
    %A = rmfield(A,{});
    
    
else
    
    % Original formatted totals
    %
    % Convert single arrays to double
    A = structfun(@double,S,'UniformOutput',0);
    
    A.names = S.names;

    
end

end


function test_case
% TEST CASE
% looks like sometimes some grid points go missing


% % % Poleward flow experiment study area
% % in_lon = [-121.5 -120.2]; % longitude boundary
% % in_lat = [34.4 35.3];   % latitude boundary
% % st = datenum(2001,8,16,17,0,0):1/24:datenum(2001,8,19,16,0,0);
% 
% % SB Channel area
%     in_lon = [-121.0 -119.0]; % longitude boundary
%     in_lat = [33.8 35.0];   % latitude boundary
% st = datenum(2004,8,16,17,0,0):1/24:datenum(2004,8,19,16,0,0);
% 
% 
% [U,V,stime,gridd,Count,AngD,Err]=get_total_data(st(1),st(end),in_lon,in_lat);
% 
% percent = total_coverage_map(U,stime,AngD,[20 160],gridd);
% 
% keyboard


% TEST HFRP FORMAT

% SB Channel area
    in_lon = [-121.0 -119.0]; % longitude boundary
    in_lat = [33.8 35.0];   % latitude boundary
st = datenum(2012,6,11,0,0,0):1/24:datenum(2012,6,11,16,0,0);


[U,V,stime,gridd,Count,AngD,Err]=get_total_data(st(1),st(end),in_lon,in_lat);

percent = total_coverage_map(U,stime,AngD,[20 160],gridd);

keyboard


end


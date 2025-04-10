function [eastkm,northkm]=rngbear2km(site_loc,central_loc,range,bearing)
% RNGBEAR2KM.M
% [eastkm,northkm]=rngbear2km(site_loc,central_loc,range,bearing)
% converts the SeaSonde radial data, given as range and bearing 
% to km East and North relative to the origin. The site location and
% central location are input as a 1x2 dim array, [longitude latitude].
% Range is in km and bearing is CCW from east.
%
% This function calls Rich Pawlowicz's dist.m which computes distances on 
% the earth using the wgs84 spheroid (using Vincenty's), not to be confused
% with the HFRProgs function of the same name.

% 22Dec97 Brian Emery 

% NEED TO MAKE SURE
% do a forward and reverse test to make sure range and bearing to km calc
% is ok - ie I need the reverse of dist.m (given range and bearing compute
% lon lat between points)
%
% see:
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=
% 15285&objectType=FILE

% reduce the confusion
site_lon=site_loc(1);
site_lat=site_loc(2);
lat_o=central_loc(2);
lon_o=central_loc(1);

% compute km east and north relative to site location
kmN_site=range.* sin(bearing.*pi./180);
kmE_site=range.* cos(bearing.*pi./180);

% compute the distance between site and origin (in two compontents)
% first the distance between the two longitudes at this lat. Note that
% the output of dist.m is always positive!
% [RANGE,AF,AR]=DIST(LAT,LONG)
[delta_E,AF,AR]=dist([lat_o lat_o],[site_lon lon_o]);

% now the distance between the two latitudes at this lon
[delta_N,AF,AR]=dist([site_lat lat_o],[lon_o lon_o]);

% compute km east and north relative to the origin
if site_lon > lon_o
	eastkm=kmE_site+(delta_E/1000);
elseif site_lon <= lon_o
	eastkm=kmE_site-(delta_E/1000);
end

if site_lat > lat_o
	northkm=kmN_site+(delta_N/1000);
elseif site_lat <= lat_o
	northkm=kmN_site-(delta_N/1000);
end

% % check that the outputs are column vectors
% if ~isempty(eastkm)
% 	eastkm=eastkm(:);
% end 
% 
% if ~isempty(northkm)
% 	northkm=northkm(:);
% end 

return

% test numbers

ptc_loc = [-120.4637 34.4543];  
cop_loc = [-119.8767 34.4100];  
central_loc = [(ptc_loc(1)+cop_loc(1))/2  (ptc_loc(2)+cop_loc(2))/2];                                  
 
	angle=181:5:361;
	range=1.5:1.5:45;
	vel=50; %this is velocity in cm/s
	current=vel*cos((angle-90)*pi/180);
	
	% make columns
	angle=angle(:);
	range=range(:);
	current=current(:);
	fudge=5*ones(size(angle));

	%use a loop to build the file: the full span of angles for each
	%range. 
	data=[];
	for i=1:length(range)
		data=[data; range(i).*ones(size(angle)), angle, current, fudge];
end


[eastkm,northkm]=rngbar2km(cop_loc,central_loc,data(:,1),data(:,2));
plot(eastkm,northkm,'r.'), hold on
[eastkm,northkm]=rngbar2km(ptc_loc,central_loc,data(:,1),data(:,2));
plot(eastkm,northkm,'b.')
title('ptc data in blue, cop data in red')
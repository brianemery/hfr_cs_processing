function LonLat = rangeBear2LonLat(RangeBear,SiteOrigin)
% RANGE BEAR TO LON LAT - mash up of other functions for a short cut
% 
% INPUT
% input bearings are usually ccwE, so assume that is the case
% 
% OUPUT
%
% EXAMPLE
%

% Copyright (C) 2017 Brian Emery
%
% Version: CamelCasE

% check for test case
if strcmp('--t',RangeBear), test_case, return, end

%   [eastkm,northkm]=rngbear2km(site_loc,central_loc,range,bearing)
%   converts the SeaSonde radial data, given as range and bearing 
%   to km East and North relative to the origin. The site location and
%   central location are input as a 1x2 dim array, [longitude latitude].
%   Range is in km and bearing is CCW from east.

[eastkm,northkm] = rngbear2km(SiteOrigin,SiteOrigin,RangeBear(:,1),RangeBear(:,2));


%   KM2LONLAT.M - Convert km to lon lat from origin using Vincenty's
%   [lon, lat] = km2lonlat_new(lon_origin,lat_origin,east,north)
%   
%   INPUT
%   decimal lon, lat of origin (scalar)
%   km east and north of origin (vectors)
%   
%   OUTPUT
%   decimal lon, lat (vectors)
 
[LonLat(:,1), LonLat(:,2)] = km2lonlat_new(SiteOrigin(1,1),SiteOrigin(1,2),eastkm,northkm);


end

function test_case
% TEST CASE

SiteOrigin = [-119.8787   34.4076];

% bearings ccwE
RangeBear(:,2) = 180:360;

RangeBear(:,1) = 30;

LonLat = rangeBear2LonLat(RangeBear,SiteOrigin);


% % now go back
% [east, north] = lonlat2km(SiteOrigin(1),SiteOrigin(2),LonLat(:,1),LonLat(:,2));
% 
%... jfc, what a hassle
for i = 1:size(LonLat,1)
    [rng(i),af(i),ar(i)] = dist([SiteOrigin(2) LonLat(i,2)],[SiteOrigin(1) LonLat(i,1)]);
    
end

% pretty close ...
keyboard


end
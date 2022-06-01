function LonLat = rangeBear2LonLat(RangeBear,SiteOrigin)
% RANGE BEAR TO LON LAT - mash up of other functions for a short cut
% LonLat = rangeBear2LonLat(RangeBear,SiteOrigin)
% 
% INPUT
% input bearings are usually ccwE, so assume that is the case
% Range in km
% 
% OUPUT
%
% EXAMPLE
%
%
% See also, compute_heading.m ...

% Copyright (C) 2017 Brian Emery

% check for test case
if strcmp('--t',RangeBear), test_case, return, end

% NOTE 
% Regarding empty radial data cases for the RNG/CS processing. It looks 
% like m_fdist.m will output all empty, if given all emtpy. But the use of
% the site lonlat here will cause an error. So, put in this trival check
if isempty(RangeBear)
    LonLat = RangeBear; return
end


%  Note that m_fdist calls for azimuth input, azimuth is compass bearing,
% eg degrees cwN, meters! too
[LonLat(:,1),LonLat(:,2),~] = m_fdist(SiteOrigin(1),SiteOrigin(2),ccwE2cwN(RangeBear(:,2)),RangeBear(:,1)*1000);

% Lon's are messed up too
LonLat(:,1) = LonLat(:,1) - 360;

return

% OLD WAY

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
RangeBear(:,2) = 160:360;

RangeBear(:,1) = 30;

LonLat = rangeBear2LonLat(RangeBear,SiteOrigin);


% % % now go back
% % [east, north] = lonlat2km(SiteOrigin(1),SiteOrigin(2),LonLat(:,1),LonLat(:,2));
% % 
% %... jfc, what a hassle
% for i = 1:size(LonLat,1)
%     [rng(i),af(i),ar(i)] = dist([SiteOrigin(2) LonLat(i,2)],[SiteOrigin(1) LonLat(i,1)]);
%     
% end
% 
% % pretty close ...
% keyboard
% 
% % 
% % bearing_f is compass bearing from 1 to 2, ie from site to grid (cwN)
% [M.Range,M.bearing_f,M.bearing_r] = m_idist(SiteOrigin(1),SiteOrigin(2),LonLat(:,1),LonLat(:,2));
% M.Bear = cwN2ccwE(M.bearing_f);


%  This is the right way (Note that azimuth is compass bearing, eg
% degrees cwN) ... less than 1e-7 difference (deg) and 1e-6 km
% [lon,lat,~] = m_fdist(SiteOrigin(1),SiteOrigin(2),ccwE2cwN(RangeBear(:,2)),RangeBear(:,1));
% 
[I.Range,I.bearing_f,I.bearing_r] = m_idist(SiteOrigin(1),SiteOrigin(2),LonLat(:,1),LonLat(:,2));
I.Bear = cwN2ccwE(I.bearing_f);
I.Range = I.Range./1000;

% these should be true:
isequal( round(I.Bear*1000)./1000, round(RangeBear(:,2).*1000)./1000 )
isequal( round(I.Range*100000)./100000, round(RangeBear(:,1).*100000)./100000 )

keyboard


% Make a map test too

mvco_map, title('click three locations')
[lon,lat] = ginput(3);


[I.Range,I.bearing_f,I.bearing_r] = m_idist(lon(1),lat(1),lon,lat);
I.Bear = cwN2ccwE(I.bearing_f);
I.Range = I.Range./1000;

% now get it back
LonLat = rangeBear2LonLat([I.Range I.Bear],[lon(1) lat(1)]);

hold on
h1 = plot(lon,lat,'o');

h2 = plot(LonLat(:,1),LonLat(:,2),'r*');

legend([h1(1) h2(1)],'original','recomputed')

keyboard

end
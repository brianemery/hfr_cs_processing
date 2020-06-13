function [Vp,bearing_r,Range,bearing_f] = uv2radcomp(site_loc,gridd,U,V)
% UV2RADCOMP - convert UV currents to radial component for HF site
% [RadComp,Bear,Range] = uv2radcomp(site_loc,gridd,U,V)
%
% Computes the component of current data in the
% direction of an HF site.
%
% INPUT
%
% OUTPUT
% Radial component of velocity (positive toward site)
% Bearing from grid to site (deg cwN)
% Range in km
% Bearing from site to grid (deg cwN)
% 
% EXAMPLE
% [Vr,Head,Range,Bear] = uv2radcomp(loc,TUV.LonLat,TUV.U,TUV.V);
% 
%
% From moordatcomp.m. See also tuv_to_radcomp.m which uses this.
 
% Copyright (C) 2014  Brian Emery
%
% version 10-Sep-2014 

% check for test case
if strcmp('--t',site_loc), test_case, return, end


% [Range,bearing_r ] = deal(NaN(size(gridd,1),1));
% 
% for i = 1:size(gridd,1)
%     
%     % bearing_r is reverse bearing in deg CWN (aka Heading?)
%     [Range(i),bearing_f(i),bearing_r(i)] = dist([site_loc(2) gridd(i,2)],[site_loc(1) gridd(i,1)]);
%     
% end

% This is essentially the same, but faster and consistent with m_fdist
% bearing_f is compass bearing from 1 to 2, ie from site to grid (cwN)
[Range,bearing_f,bearing_r] = m_idist(site_loc(1),site_loc(2),gridd(:,1),gridd(:,2));

% convert to km
Range = Range/1000;

% convert to radians
br_rad = radians(bearing_r); 


% Expand bearins into matrix if necessary
if size(U,2) > 1
    br_rad = repmat(br_rad,1,size(U,2));
end


% Use the bearing to compute the component of velocity radial to the codar
% note that positive v_parallel is in the direction of the radar
% This is the vector dot product, see arof.m notes in notebook for math.
 Vp = (U.*sin(br_rad))+(V.*cos(br_rad));
 


end




function test_case

% SIMPLE TEST
site = [-120.4716 34.4483];

u = 50; v = 20;

gridd = [-120.6543   34.2713];

[Vp,Brg] = uv2radcomp(site,gridd,u,v);

keyboard


sbc
hold on
site
plot(site(1),site(2),'r*')
plot(gridd(1),gridd(2),'r*')




return
% CONVERT THIS INTO A TEST

R.SiteOrigin = [-120.4716 34.4483];

% load totals for whole day
load /projects/error_covariance/data/tot_2013_05_27.mat

% generate all the range cell locations
rng = unique(R.RangeBearHead(:,1));

% brg ccwe
bear = [100:350]'; rng = rng*ones(size(bear));

[ekm,nkm] = rngbear2km(R.SiteOrigin,R.SiteOrigin,rng,bear);

[lon,lat] =km2lonlat(R.SiteOrigin(1),R.SiteOrigin(2),ekm,nkm);


% interp UV data to the range cell locations
[u,v,er,TUV] = interp_to_codar(lon,lat,R.TimeStamp,TUV);


% compute radial component
RadComp = uv2radcomp(R.SiteOrigin,[lon lat],u,v);
% [~,RadComp] = moordatcomp(R.SiteOrigin,[lon lat],u,v);



end
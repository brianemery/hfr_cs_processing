function [lon,lat] = km2lonlat_new(lon_o,lat_o,x,y) 
% KM2LONLAT.M - Convert km to lon lat from origin using Vincenty's
% [lon, lat] = km2lonlat_new(lon_origin,lat_origin,east,north)
% 
% INPUT
% decimal lon, lat of origin (scalar)
% km east and north of origin (vectors)
% 
% OUTPUT
% decimal lon, lat (vectors)
%
% Basically an attempt to reproduce the output of Codar's total gridding
% program. Difference with theirs is ~25 cm over 30 km (mostly in x 
% direction) at high latitudes (~60 N in test case).
%
% uses the M_Map tools m_idist.m and m_fdist.m, which use the WGS84 earth
% model and Vincenty's algorithm.
%
% See also test_km2lonlat_method, lonlat2km.m

% Copyright (C) 2010 Brian M. Emery
% Version 1.0, 1 Oct 2010
%              5 Oct 2010 debugging logical indexing in subfunction

% Optional test case
if strcmp(lon_o,'--t')
    test_case, [lon,lat] = deal(NaN); return
end
  
% reshape, remove nans
[x,y,r,c,i] = reform_inputs(x,y);

% Convert xy to meters
x = x.*1000;
y = y.*1000;

% Expand origin input into matrix 
lon_o = lon_o.*ones(size(x));
lat_o = lat_o.*ones(size(x));

% define the positive direction bearings
xbrg = 90 .* ones(size(x));
ybrg = zeros(size(x));


% use y to get lat for x calc
[lon2,lat2] = m_fdist(lon_o,lat_o,ybrg,y);

% apply x to get lon lat 
[lon,lat] = m_fdist(lon2,lat2,xbrg,x);

% x = 0 case
lon(x == 0) = lon_o(x == 0);

% Convert longitude convention
lon(lon > 180) = lon(lon > 180) - 360;

% reinsert nan and reshape
[lon,lat] = reform_outputs(lon,lat,r,c,i); 

end
%% ----------------------------------------------------------------
function [x,y,r,c,i] = reform_inputs(x,y)
% REFORMAT INPUTS
% reshape to columns, remove nans

% reshape to columns, keep r,c
[r,c] = size(x+y);
x = x(:);
y = y(:);

% remove nans, reinsert them in the output
i = ~isnan(x+y);
x = x(i);
y = y(i);

end
%% ----------------------------------------------------------------
function [ln,lt] = reform_outputs(lon,lat,r,c,i)
% REFORM OUTPUTS

% % Potential matlab bug where this sometimes generates zeros
% % reinsert Nan's
% lon(~i) = NaN;
% lat(~i) = NaN;

% work around:
[ln,lt] = deal(NaN(r*c,1));
ln(i) = lon;
lt(i) = lat;

% reshape outputs
ln = reshape(ln,r,c);
lt = reshape(lt,r,c);

end


%% ----------------------------------------------------------------
function test_case
% TEST CASE for the development of this lonlat2km version
%
% Construct the lon lat grid using m_map tools, then compute the km
% equivalents.

% ----------------------------------------------------------------
% BASIC TEST CASE
% ----------------------------------------------------------------

x = NaN(6,111);
x(:,1) = 1:6;
y = x;

% run the function
[lon,lat] = km2lonlat(-140,60,x,y);
keyboard


% ----------------------------------------------------------------
%  COMPARE VS COS TOTAL FILE
%------------------------------------------------------------------

% LonLat_grid uses m_map tools to construct a grid. This results in a grid
% that is very similar to eg a grid found in Codar's total files. An
% example can be found here:
load /projects/drifter_simulation/pws/data/totals/tot_pws.mat

% the above contains lon lats from the Prince William Sound, which provides a pretty
% robust test, and XY from the Codar total vector gridding program.

% the origin is here:
ii = find(TUV.OtherSpatialVars.X ==0 & TUV.OtherSpatialVars.Y ==0);

[lon,lat] = km2lonlat(TUV.LonLat(ii,1),TUV.LonLat(ii,2),TUV.OtherSpatialVars.X,TUV.OtherSpatialVars.Y);


% CHECK PLOTS

% grid comparison
compare_grids(lon,lat,TUV.LonLat(:,1),TUV.LonLat(:,2))

% compute and plot the differences
plot_ll_diffs(TUV.LonLat(:,1),TUV.LonLat(:,2),lon,lat);
title('COS vs km2lonlat total differences (meters)')




% ----------------------------------------------------------------
% COMPARE WITH lonlat2km_dev and back again
%------------------------------------------------------------------

% Get XY
[x,y] = lonlat2km(TUV.LonLat(ii,1),TUV.LonLat(ii,2),TUV.LonLat(:,1),TUV.LonLat(:,2));

% GET LONLAT BACK
[lon,lat] = km2lonlat(TUV.LonLat(ii,1),TUV.LonLat(ii,2),x,y);

% CHECK PLOTS

% grid comparison
compare_grids(lon,lat,TUV.LonLat(:,1),TUV.LonLat(:,2))

% compute and plot the differences
plot_ll_diffs(TUV.LonLat(:,1),TUV.LonLat(:,2),lon,lat);
title('lonlat2km then km2lonlat total differences (meters)')



% % ----------------------------------------------------------------
% % COMPARE WITH LonLat_grid.m
% %------------------------------------------------------------------
% % origin is a bit of a problem ...
% % these from PWS grid
% LowerLeft  = [-147.5181635   60.3050916];
% UpperRight = [-146.4599690   60.7720839];
% 
% %
% [ Lon, Lat ] = LonLat_grid( LowerLeft, UpperRight, 2);
% 
% keyboard
% x = -28:2:30;
% y = -28:2:24;
% 
% [lon,lat] = km2lonlat(lon_o,lat_o,x,y) 


end
%% ----------------------------------------------------------------
%  TEST CASE PLOTTING
function compare_grids(x,y,a,b)

figure
plot(x,y,'r.')
hold on
plot(a,b,'bo')
title('Grid Comparison')

end
%% ----------------------------------------------------------------
function plot_ll_diffs(lon1,lat1,lon2,lat2)

[s,a12,a21] = m_idist(lon1,lat1,lon2,lat2);

figure
cdot2d(lon1,lat1,s);


end
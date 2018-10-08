 function [Axis_factor, Scale_factor] = mercat(lons,lats)

%MERCAT  mercator scaling factors for axes 'AspectRatio' parameter.
%  This function will calculate scale factors to be used in the axes
%  command option 'AspectRatio' to produce mercator projection plots.
%  The latitude used to scale the longitudinal distance is simply 
%  the midpoint of the 2 input latitudes.
%
%  The format is as follows:
%  [Axis_factor, Scale_factor] = mercat(lons,lats)
%
%  input:
%         lons is a 2 element vector - [minimum lon, maximum lon]
%	  lats is a 2 element vector - [minimum lat, maximum lat]
%
%  output:
%         Axis_factor  - scalar
%         Scale_factor - scalar
%
%  usage example:
%         axes(... ,'AspectRatio',[Axis_factor,Scale_factor], ...)
%
%  See the "axes" command discussion in the MATLAB Reference Guide
%  for more information on AspectRatio.
%

% 	Mike Cook - NPS Oceanography Dept. -    Oct 93
%   Added error statements                  Mar 94

 if nargin ~= 2
    error(' Must supply 2 element longitude and latitude vectors ... type help mercat')
 end
 if length(lons) ~= 2  |  length(lats) ~= 2
    error(' Must supply 2 element longitude and latitude vectors ... type help mercat')
 end

 lat_mid = (lats(2) + lats(1)) / 2;
 Scale_factor = cos((lat_mid*pi)/180);

 Axis_factor = ( (lons(2) - lons(1)) / (lats(2) - lats(1)) ) * Scale_factor;


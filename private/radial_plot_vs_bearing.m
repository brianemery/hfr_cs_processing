function h = radial_plot_vs_bearing(R,rcell,tdx,clr)
% RADIAL PLOT VS BEARING - plot radial vs bearing
% h = radial_plot_vs_bearing(R,rcell,tdx,clr)
% 
% INPUT:
% RADIAL struct
% range cell index
% time index
% line style, color (eg '-r.')

% Copyright (C) 2010 Brian M. Emery
% 22 June 2010

if nargin < 4, clr = 'bo'; end

% Get unique ranges
rngs = unique(R.RangeBearHead(:,1));

% find the one of interest
j = find(R.RangeBearHead(:,1) == rngs(rcell));

% plot it
h = plot(R.RangeBearHead(j,2),R.RadComp(j,tdx),clr);


end

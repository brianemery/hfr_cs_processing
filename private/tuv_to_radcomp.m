function [Vr,Bear,Range] = tuv_to_radcomp(TUV,loc,n)
% TUV TO RADCOMP - radial component of totals, sorted by range?
% Vr = tuv_to_radcomp(TUV,loc,n)
% 
% This is different than other similar mfiles because it takes the RADIAL
% field as unknown and creates it from the TUV data along with info about
% the bins to make
%
% ... Then make it a tool ...  we're not really
% making radials here , just the radially binned data ... 
%
% SEE ALSO
% makeRadialsFromTotals.m, uv2radcomp (used here)

% TO DO?
% store in RADIAL format? Keep LonLats, RangeBearHead, UV??


% ***  JUST GRAB ONE RANGE CELL FOR NOW ***
%
% VERY SLOW ... screen by range first?


% Get radial component of velocity ... and other stuff easily computed at
% the same time
[Vr,Head,Range,Bear] = uv2radcomp(loc,TUV.LonLat,TUV.U,TUV.V);




i = find( Range < 21 & Range > 19.5);

keyboard

% sbc
% plot(lon(i),lat(i),'b.')
% 




% Vp(Vp==0) = NaN;
% 
% plot(br,Vp*100,'.')
% 
% xlabel('Bearing (^o CWN)')
% ylabel('Radial Velocity (cm s^-^1)')
% 
% keyboard
% 
% publicationStandards
% saveas(gcf,'~/Dropbox/roms_radial_component_figure_3.png')
% 
% 



end


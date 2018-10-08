function maskIdx=sbgrid_mask(gridd)
% SBGRID_MASK.M
% maskIdx=sbgrid_mask(gridd);
% Creates a variable called maskIdx which contains the indicies of the
% grid points from sbgrid.m which are on or near land. 
%
% Usage:
%
% % blank the grid points in the total vector matrix by setting them
% % to NaN's
% maskIdx=sbgrid_mask(gridd)
% totalUVm(maskIdx,:)=NaN;
%
% Or,
% U(maskIdx,:)=NaN;
% V(maskIdx,:)=NaN;
%
% However, this is now called by sbgrid, so these grid points may not be in
% the final grid.

% Brian Emery July 2008 major updates and improvements

if length(nargout)==0
    disp('Must Specify an output for sbgrid_mask.m. Suggest:')
    disp('maskIdx=sbgrid_mask(gridd);')
    keyboard
end

% initialize 
mask(1:length(gridd),1)=1;

% Use Coastline file data to make the polygon. 'k' is the indecies of each
% of the polygons ... 
coastFile='/Data/Mapping/COAST4_124_29.mat';
eval('load(coastFile);','error([coastFile '' not found''])')

% trial and error sets n as the indicies of mainland, SCI, SRI, SMI and SNI
for n=[1 2 5 10 12]
    % this gives a logical 'in', ==1 if the point is in the polygon 
    in=inpolygon(gridd(:,1),gridd(:,2),ncst(k(n):k(n+1),1),ncst(k(n):k(n+1),2));
    mask(in==1)=NaN;
end

maskIdx=find(isnan(mask));

% keyboard
% % check plot
% sbc, plotgrid
% plot(gridd(maskIdx,1),gridd(maskIdx,2),'c.')
% 
% % also eliminate some points in the water but close to land
%  maskIdx=unique([maskIdx(:)' 22:33 38:47 68:76 84:94 109:112 119:123 132:141 153:158 170 178:186 201:204 ...
% 224:231 250 271:273 ...
% 1125:1128 1165:1175 1209:1222 1232:1235 1254:1269 1279:1316 1326:1363 1372:1410 ...
% 1419:1457 1464:1504  1511:1551]);

% % keep unique indicies
% maskIdx=unique(maskIdx);
% 
end

%% ---------------------------------------------------------------------
% Useful code not used here ...
function TUV = screen_and_clean(TUV)

% a = axis to get box
a = [-120.3398 -119.2032   34.0070   34.4773];


% this gives a logical 'in', ==1 if the point is in the polygon 
in=inpolygon(TUV.LonLat(:,1),TUV.LonLat(:,2),a([1 2 2 1 1]),a([3 3 4 4 3]));

keyboard
    
TUV = subsref_struct(TUV,find(in),size(TUV.LonLat,1),1);


end
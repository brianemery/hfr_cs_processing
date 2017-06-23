function [Lon,Lat,Z,address] = my_griddata(lon,lat,z)
% MY_GRIDDATA
% [LON,LAT,Z,address] = my_griddata(lon,lat,z)
% 
% take z,lat,lon column vectors and put them into lat lon gridded matricies
% similar to griddata, but no interpolation. See also address_matrix.m 
%
% z can be a matrix, with # rows = length(lon) and columns for each time
% (for example). Output Z will be 3-d
%
% The vector inputs must represent a regular grid for this to work
% correctly.
%
% outputs optional address index which could be used like this:
%
% [Lon,Lat,U,address] = my_griddata(lon,lat,u);
% V = NaN(size(U)); 
% V(address) = u;

% a breathtaking work of staggering genius from stuffForTonyaGnomeFormat2.m
% Brian Emery ca 2007

% orient inputs
lon = lon(:);
lat = lat(:);
if size(z,1) ~= length(lon)
    disp('input size error'), keyboard
end

% put data into matrix
[latQ,ix,rowIdx]=unique(lat); clear ix
[lonQ,iy,colIdx]=unique(lon); clear iy
[Lon,Lat]=meshgrid(lonQ,latQ);

% convert row,col indecies into matrix addresses
address=( size(Lat,1).*(colIdx -1) )+ rowIdx ;

% create the Z matrix
% Z=NaN(size(Lat));
Z = NaN(size(Lat,1),size(Lat,2),size(z,2));

% fill in data. Not sure how to do this w/out loop
for t = 1:size(z,2)
    zgrd = NaN(size(Lat));
    zgrd(address) =z(:,t);
    Z(:,:,t)=zgrd;
end

% % VERIFY THIS IS CORRECT
% cdot2d(lon,lat,z)
% hold on
% figure
% cdot2d(Lon,Lat,Z)
% for j = 1:10, figure(1), pause,figure(2), pause, end

% ALTERNATE METHOD
% keyboard
% 
% Z=NaN(size(Lat));
% 
% for i = 1:length(z),    Z(Lat(:,1) == lat(i), Lon(1,:) == lon(i)) = z(i); end
% 

end

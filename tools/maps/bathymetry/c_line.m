 function Cline=c_line(lon,lat,data,level)

% C_LINE.M
% [Cline]=c_line(x,y,z,level)
% Outputs 'Cline', the x and y positions
% of a given contour line given in 'level' so that the data 
% can be used with 'plot' to plot the contour. The x,y 
% inputs must be vectors describing the positions of the
% data given in the matrix z, as in the 'contour' function.
%
% example:
% load c:\data\bathymetry\sbcsmb.mat
% [d200]=c_line(lon,lat,depth,200);
% plot(d200(:,1),d200(:,2))
%
% See code at the end of sbgrid2.m for creating contourable
% matricies on the HF grid.

% try [d700]=c_line(lon',lat,depth,700); with sbcsmb_lo.mat

% Brian Emery 16Apr99

% compute the contour matrix 
C = contourc(lon,lat,data,[level level]);

% define m, the index of the contour level in C(see help contourc)
m=find(C(1,:)==level);

Cline=C;
Cline(:,m)=NaN;
Cline=Cline';

return

% this quickly makes a bathy data set
%for i=[700:100:1000 1200:200:3600];
for i=[325:50:500];
x=num2str(i);
eval(['[d' x ']=c_line(lon,lat,depth'',' x ');'])
eval(['plot(d' x '(:,1),d' x '(:,2),''b'')']), hold on

end
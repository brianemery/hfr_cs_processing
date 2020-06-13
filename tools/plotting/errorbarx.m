function h = errorbarx(x, y, l,u,symbol)
% ERRORBARX.M
% Same as Errorbar.m except the bar is drawn in x instead of y, 
% and the errorbars only overplot onto an existing plot.
%
% Help for regular errorbar:
% 
%     ERRORBAR(X,Y,L,U) plots the graph of vector X vs. vector Y with
%     error bars specified by the vectors L and U.  L and U contain the
%     lower and upper error ranges for each point in Y.  Each error bar
%     is L(i) + U(i) long and is drawn a distance of U(i) above and L(i)
%     below the points in (X,Y).  The vectors X,Y,L and U must all be
%     the same length.  If X,Y,L and U are matrices then each column
%     produces a separate line.
%  
%     ERRORBAR(X,Y,E) or ERRORBAR(Y,E) plots Y with error bars [Y-E Y+E].
%     ERRORBAR(...,'LineSpec') uses the color and linestyle specified by
%     the string 'LineSpec'.  The color is applied to the data line and
%     error bars while the linestyle and marker are applied to the data 
%     line only.  See PLOT for possibilities.
%  
%     ERRORBAR(AX,...) plots into AX instead of GCA.
%  
%     H = ERRORBAR(...) returns a vector of errorbarseries handles in H.
%  
%     For example,
%        x = 1:10;
%        y = sin(x);
%        e = std(y)*ones(size(x));
%        errorbar(x,y,e)
%     draws symmetric error bars of unit standard deviation.
% 
%     Reference page in Help browser
%        doc errorbar
% 

% Changes made 6jun96 - Brian Emery 

if min(size(x))==1,
  npt = length(x);
  x = x(:);
  y = y(:);
  if nargin == 3,  
			l = l(:);
		elseif nargin == 4 | nargin == 5
			l = l(:);
			u = u(:);
		end
else
  [npt,n] = size(x);
end

if nargin == 3,  
	u = l;
end

if nargin ~= 5
			symbol = '-';
end

if nargin == 2
  l = y;
		u = y;
  y = x;
  [m,n] = size(y);
  x(:) = [1:npt]'*ones(1,n);;
end
if isstr(x) | isstr(y) | isstr(l) | isstr(u)
	error('Arguments must be numeric.')
end

if any(size(x)~=size(y)) | any(size(x)~=size(l)) |  any(size(x)~=size(u)),
  error('The sizes of X, Y, L and U must be the same.');
end

% % This is how matlab did it, which is fine usually
% tee = (max(y(:))-min(y(:)))/100;  % make tee .02 y-distance for error bars
% ... but for layering you get different size tee's. So base it on the axis
% instead 
a = axis; tee = (a(4)-a(3))/100;

yl = y - tee;
yr = y + tee;
n = size(x,2);

% Plot graph and bars
cay = newplot;
next = lower(get(cay,'NextPlot'));
% build up nan-separated vector for bars
xb = [];
yb = [];
nnan = nan*ones(1,n);
for i = 1:npt
    xtop = x(i,:) + l(i,:);
    xbot = x(i,:) - u(i,:);
    yb = [yb; y(i,:); y(i,:) ; nnan; yl(i,:);yr(i,:);nnan ;yl(i,:);yr(i,:);nnan];
    xb = [xb; xtop;xbot;nnan;xtop;xtop;nnan;xbot;xbot;nnan];
end

h = plot(xb,yb,symbol);


end

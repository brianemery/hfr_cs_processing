function [ax,p] = cdot2d(x,y,cval,dotsize,MapColor,crange)
%CDOT2D  Draws a 2-D pseudocolor plot w/data represented as points.
%
%  CDOT2D(X,Y,CVAL) creates a 2-D plot with the elements of vector CVAL
%  represented as pseudocolored dots plotted at location X, Y.
%  The values of the elements of CVAL are mapped to one of Matlab's
%  pseudocolor tables.  A colorbar is also plotted.
% 
%  NOTE:  X, Y, and CVAL must all be row or column vectors of the same
%  length!  Matrices will not work properly.
%
%  You MUST supply the first 3 arguments, the following 4 are
%  optional and if not supplied will be set to default values.
%
%  CDOT2D(X,Y,CVAL,DOTSIZE) will plot the CVAL data as dots of size DOTSIZE.
%  The default DOTSIZE is set to 20.
%
%  CDOT2D(X,Y,CVAL,DOTSIZE,CMAP) will plot the data using the color map
%  CMAP.  The default color map is JET.  Type 'help graph3d' for a list
%  of Matlab's color maps.  CMAP can be a string, or an n by 3 matrix
%  defining the colormap.
%
%  CDOT2D(X,Y,CVAL,DOTSIZE,CMAP,CRANGE) where CRANGE is a 2 element vector
%  [crange_min, crange_max], limits the range of the pseudocolored plot from
%  crange_min to crange_max.  CVAL data outside this range are not plotted.
%  The default is to use the min and max values of CVAL for CRANGE.
%
%	Example:
%
%	CDOT2D(randn(30,1),randn(30,1),randn(30,1),25,'copper',[-2,2]);
%	would plot 30 random points using dots with a size of 25,
%       and the copper color map with a data range from -2 to 2.
%
%  Help graph3d will give you a listing of matlab's built in color maps.
%
% BE MODS:
% [ax,p] = cdot2d ... outputs colorbar axis handle and handles for the many 
%          dots.
%
% SEE ALSO: colordot.m


%	Mike Cook - NPS Oceanography Dept., JULY 94
%
%   Updates by Brian Emery, Nov 2013

%%%%%%  Default assignments and error checking section %%%%%%%
dfltdotsize =20;		% Default dot size
dfltMapColor = 'jet';		% Default colormap

if nargin < 3
   error(' You *MUST* supply at least 3 arguements');
end
if nargin < 4
   dotsize = dfltdotsize;
end

if nargin < 5
   MapColor = dfltMapColor;
end

if nargin < 6	% If user doesn't supply max/min range - use max/min of data.
   crange(1) = min(min(cval));
   crange(2) = max(max(cval));
else
   if length(crange) ~= 2
     error(' You *MUST* supply a 2 element vector for the data range')
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% enable input of map as colormap array
if ischar(MapColor)
    colormap(MapColor);
    map=colormap;			% Get color scale - 64 rgb triplets.
else
    map = MapColor;
end
    

x = x(:);
y = y(:);
cval = cval(:);

row= length(cval);

b1 = (size(map,1)-1) / (crange(2)-crange(1));
b0 = ( -b1 * crange(1) ) + 1;
Clrs = round( (cval .* b1) + b0);

% % % Draw points with the value of each point represented by a color.
% % mainAX = axes('Box','on','TickDir','out','TickLength',[0.017,0.025], ...
% %              'FontWeight','bold');

p = NaN(row,1);

% BE: could vectorize this, just not at the moment
              
for i = 1:row
    if Clrs(i) >= 1  &&  Clrs(i) <= 64
%         % Code to plot dots, or ...
%         p(i) = line('Xdata',x(i),'Ydata',y(i),'Marker','.', ...
%                         'MarkerSize',dotsize,'Color',map(Clrs(i),1:3) );

        % ... code to plot squares.
        line('Xdata',x(i),'Ydata',y(i),'Marker','s', ...
            'MarkerSize',dotsize,'Color',map(Clrs(i),1:3), ...
            'MarkerFaceColor', map(Clrs(i),1:3), ...
            'MarkerEdgeColor','none');

        %   else
        %     line('Xdata',x(i),'Ydata',y(i),'Marker','.', ...
        %        'MarkerSize',dotsize,'Color','w' );
    end
end


% % BE DISABLE THIS STUFF
% % might create side effects, need to check, but I want this to be a bit
% % more general purpose going forward?
% 
% mainAX = gca;
% set(mainAX,'Box','on','TickDir','out','TickLength',[0.017,0.025], ...
%            'FontWeight','bold');
% 
%        keyboard
%        
% % resize surface object axis - preserve aspectratio if set.  If aspectratio
% % is not set, the axis and scale factors are set to [NaN, NaN].
% aspectfactors = get(mainAX,'DataAspectRatio');
% 
% 
% % return % removes colorbar
% 
% % % Draw a colorbar which shows the scale of the points.
% % % resize main axis so colorbar fits vertically to the left of the main axis.
% % set(mainAX,'units','normal');
% % opos=get(mainAX,'pos');
% % dx = opos(3)/(.775/.15);
% % set(mainAX,'Position',[opos(1)+dx, opos(2), opos(3)-dx, opos(4)], ...
% %            'DataAspectRatio',aspectfactors);
% % %%set(mainAX,'Position',[opos(1)+0.15, opos(2), opos(3)-0.15, opos(4)]);
% % 
% % % create new axis for colorbar
% % dx = opos(3)/(.775/.01);
% % dx2 = opos(3)/(.775/.05);
% % ax = axes('Position',[opos(1)-dx, opos(2), dx2/2, opos(4)],'units','normal');
% % %%ax = axes('Position',[opos(1)-0.01, opos(2), 0.05, opos(4)],'units','normal');
% 
% 
% 
% % BE HACK
% % put the colorbar on the right axis
% 
% % "Setting the DataAspectRatio will disable the stretch-to-fill behavior"
% set(mainAX,'DataAspectRatio',aspectfactors);

% % get position of main axis in normalized units
% set(mainAX,'units','normal');
% opos = get(mainAX,'pos');
% 
% % get width of colorbar based on this
% dx = opos(3)/30;
% 
% % create new axis for colorbar on right axis
% ax = axes('Position',[opos(1)+opos(3) opos(2) dx opos(4)],'units','normalized');
% 
% % Create Colorbar
% y = linspace(crange(1),crange(2),100);
% 
% ColorRange=[y', y'];
% 
% p = pcolor([1,2], ColorRange(:,1), ColorRange);
% 
% set(p,'FaceColor','interp','EdgeColor','none');

% why not use this?
%ax = colorbar('peer',mainAX,'Location','EastOutside');
mainAX = gca;
ax = colorbar_my(mainAX);

caxis(crange);

% set(ax,'Xlim',[1,2], ...
%        'Ylim',crange,...
%        'FontWeight','bold', ...
%        'TickLength',[0,0],...
%        'XTickLabel',[' '],...
%        'YAxisLocation','right', ...
%        'Box','off'); %,'YTickLabel',[' ']) %rm's colorbar labs
   
% make surface object axis default
axes(mainAX);


end

function Clrs = get_cmap_index(cval,map,crange)
% GET CMAP INDEX
% 
% weird code from cdot2d.m

% this used eslewhere, incorp at some point


cval = cval(:);

b1 = (size(map,1)-1) / (crange(2)-crange(1));

b0 = ( -b1 * crange(1) ) + 1;

Clrs = round( (cval .* b1) + b0);



end
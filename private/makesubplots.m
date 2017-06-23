function [haxes]=makesubplots(nrows,ncols,sh,sv,varargin)
% MAKE SUBPLOTS -  closely spaced subplots
% 
% [HAXES] = makesubplots(NROWS,NCOLS,sh,sv,extras)
%
% Sets up an  NROWS by NCOLS array of closely-spaced subplots in the
% current window.  Returns a length NCOLS*NROWS vector HAXES of handles 
% to the axes in column-wise order. SH and SV are the normalized 
% spacing between subplots.
% NOTE: sh=sv=0.01 for axes without labels
%       sh=sv=0.1  for axes with labels
%
% EXAMPLE:
% Make 7 plots on one figure (arranged in 7 rows and 1 column wide):
% 
% % Create the subplot axes
% % and add whatever you want to plot to it:
% hx = makesubplots(7,1,.05,.05);
% h = plot(hx(1), doy, WS, ln);
% 
% % Add letter (see subplot_add_letters.m! run it outside the loop)
% subplot_add_letters(hx)
%
% % Instead of axes.m, use:
% set(gcf,'CurrentAxes',hx(2))
%
%
% SUPPORTS FRACTIONAL SUBPLOT SIZES
% EXAMPLE:
% % Fractions total to 1
% hx = makesubplots(9,2,.04,.03, [1/9*ones(1,6) .5/9 1/9 1.5/9]);




% Bill Shaw's code
% Modified 8 Feb 2010 by Brian Emery 
%   - added fraction plot options

% TO DO:
% varargin bw compatible, expand to width subplots (?), check inputs ...
% spacing between plots should be same in fractional case (right now its 
% a function of the plot height which changes ..)
%
% Need to be able to set fractional width also
%
% automate adding letters to the plots, use this:
% char(97:97+25)

extras = {};

% define borders
bh = 0.075; %0.125; %0.075;
bv=0.075;

% if length(varargin)
%     if ~isempty(varargin{1})
%         bh = varargin{1};
%     end
% end
% if length(varargin)>2
%     if ~isempty(varargin{1})
%         bv = varargin{2};
%     end
% end
% if length(varargin)>2
%     extras = varargin(3:end);
% end



% define borders and cells and spacing
wcell=(1-2*bh)/ncols;

if ~isempty(varargin)
    hfrac = varargin{1}; % height fraction totaling 1
    hcell = (1-2*bv)*hfrac;
else
    hcell=(1-2*bv)./(nrows*ones(1,nrows));
end

% create matrix of handles
haxes = zeros(nrows,ncols);

% % hack to allow variable horiz spacing
% if length(sh) == 1, sh = sh*ones(1,ncols); end

for j=1:ncols;
    
    
	for i=1:nrows
		
        % compute address
        k = (j-1)*nrows+i;
        
        % define the axes position (left bottom width height)
        % and Add in spacing between plots
		lbwh = [ (j-1)*wcell + bh + sh*wcell, ...
                 (1-bv)-(sum(hcell(1:i))) + sv*hcell(i), ...
                            wcell - 2*sh*wcell, ...
                            hcell(i) - 2*sv*hcell(i)];
		
		haxes(i,j)=axes('Units','Normalized','Position',lbwh);
        
        % % add letter to plot
        %text(.015,0.92,[char(96+k) ')'],'units','normalized')
	end
end

set(haxes,'Box','On','NextPLot','Add',extras{:})

end
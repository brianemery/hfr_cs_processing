function publicationStandards(lbwh,ismap)
% PUBLICATION STANDARDS.M
% publicationStandards(lbwh)
% A comprehensive list of things to change 
% to make pub quality figs.
%
% Optionally input figure position (inches) 
%
% EXAMPLE
% % lets say you have a map ...
% lbwh = [0.5871    2.7218    9.8198    7.9919];
% ismap = 1; 
% publicationStandards(lbwh,ismap)

% Copyright (C) 2000-2010 Brian Emery

% NOTES
% L&O standards:
% See: http://www.aslo.org/lomethods/instructions/manuscript.html#figures
% - Tiff (compressed or not)
% - Color figures must be set to CMYK mode, with a minimum resolution of 350 dpi.
% - Line art generally should be set to a higher resolution, e.g. 1200 dpi. 
% - Grayscale images should be set to 350 dpi or higher.
% - All text and numerals on the figures are in the Times New Roman font
%   and are at least 8 points (0.1 inch, 2.8 mm) high after reduction to
%   the size that they will appear in the journal.
% - Figure legends (one paragraph per figure) explain all panels (A, B, ...).
%   Symbols used in the figure (e.g., circles, squares, ...) are explained in
%   a symbol table on the figure itself, and not in the legend. The symbol table
%    should not be surrounded by a box
% - PC or Macintosh versions of Adobe Postscript fonts should be used to label figures.
%    Do not use TrueType fonts or system "bitmap" fonts.
% - "Thus, the maximum size for a figure is 18.4 x 23.2 cm = 7-1/4 x 9-1/8 
%    inches = 2535 x 3196 pixels at 350 dpi resolution. However, some room must
%    be allowed for the figure legend, which will be displayed under the figure.
%
% AGU E-LIGHTNIING
% "1920px wide and 1080px high" png



% set(hfig,'DefaultAxesFontSize',14)
% set(hfig,'DefaultAxesLineWidth',1.5) 

if nargin < 2, ismap = false; end

fntName = 'Times New Roman';


% keep handle of current axes
htop = gca;

% Get all axes handles
haxes = findobj(get(gcf,'Children'),'type','axes'); 

% % get rid of colorbars and legends
% haxes = haxes(~strcmp('Colorbar',(get(haxes,'Tag'))));
% haxes = haxes(~strcmp('legend',(get(haxes,'Tag'))));

% Apply to whole figure:
% set(findobj('LineWidth',.5) ,'LineWidth',1.5)
set( findobj('FontSize',10)  ,'FontSize',14)

if ~ismap, 
    set(findobj(get(gcf,'Children'),'LineWidth',.5) ,'LineWidth',1.5),
else
    % gets site names
    set(findobj(get(gcf,'Children'),'FontSize',9), 'FontSize',14)
end




% Loop over each axes and change settings
for j = 1:length(haxes)
    
    % if j==length(haxes), keyboard, end
    
    % Removing this keeps all the layers in the right order
    % axes(haxes(j))
    % set(gcf,'CurrentAxes',haxes(j))
    
    % set(    findobj('LineWidth',.5) ,'LineWidth',1.5)
    % set(    findobj('FontSize',10)  ,'FontSize',14)

        
    set(    get(haxes(j),'XLabel')  ,'FontSize',14, 'FontName',fntName)   %,'FontWeight','bold')
    set(    get(haxes(j),'YLabel')  ,'FontSize',14, 'FontName',fntName)   %,'FontWeight','bold')
    set(    get(haxes(j),'title')   ,'FontSize',14, 'FontName',fntName)    %,'FontWeight','bold')
    set(    haxes(j)                ,'FontSize',14, 'FontName',fntName)
    set(    haxes(j)                ,'LineWidth',1) % 1.5


    
end


if nargin < 1
    % these are kind of arbitrary, but easy to change
    lbwh = [0.5871    2.7218    9.8198    7.9919];
end


set(gcf,'PaperOrientation','Portrait')
set(gcf,'Units','Inches')         
set(gcf,'Position',lbwh)
set(gcf,'PaperPosition',lbwh)

% Weird bug in axes.m
% axes(htop);
% Try this:
set(gcf,'CurrentAxes',htop);

set(findall(gcf,'Type','text'),'FontName',fntName)



end

function publicationStandards_other(lbwh)
% PUBLICATION STANDARDS.M
% publicationStandards(lbwh)
% A comprehensive list of things to change 
% to make pub quality figs.
%
% Optionally input figure position (inches) 

% Copyright (C) 2000-2010 Brian Emery

% NOTES
% L&O standards:
% See: http://www.aslo.org/lomethods/instructions/manuscript.html#figures
% - Tiff (compressed or not)
% - Color figures must be set to CMYK mode, with a minimum resolution of 350 dpi.
% - Line art generally should be set to a higher resolution, e.g. 1200 dpi. 
% - Grayscale images should be set to 350 dpi or higher.
% - All text and numerals on the figures are in the Times New Roman font
%   and are at least 8 points (0.1 inch, 2.8 mm) high after reduction to
%   the size that they will appear in the journal.
% - Figure legends (one paragraph per figure) explain all panels (A, B, ...).
%   Symbols used in the figure (e.g., circles, squares, ...) are explained in
%   a symbol table on the figure itself, and not in the legend. The symbol table
%    should not be surrounded by a box
% - PC or Macintosh versions of Adobe Postscript fonts should be used to label figures.
%    Do not use TrueType fonts or system "bitmap" fonts.
% - "Thus, the maximum size for a figure is 18.4 x 23.2 cm = 7-1/4 x 9-1/8 
%    inches = 2535 x 3196 pixels at 350 dpi resolution. However, some room must
%    be allowed for the figure legend, which will be displayed under the figure.


% set(hfig,'DefaultAxesFontSize',14)
% set(hfig,'DefaultAxesLineWidth',1.5) 

fntName = 'TimesNewRoman';


% keep handle of current axes
htop = gca;

% Get all axes handles
haxes = findobj(get(gcf,'Children'),'type','axes'); 

% get rid of colorbars and legends
haxes = haxes(~strcmp('Colorbar',(get(haxes,'Tag'))));
haxes = haxes(~strcmp('legend',(get(haxes,'Tag'))));

% Loop over each axes and change settings
for j=1:length(haxes)
    
    % Removing this keeps all the layers in the right order
    % axes(haxes(j))
    
    set(    findobj('LineWidth',.5) ,'LineWidth',1.5)
    set(    findobj('FontSize',10)  ,'FontSize',14, 'FontName',fntName)
    set(    get(haxes(j),'XLabel')  ,'FontSize',14, 'FontName',fntName)   %,'FontWeight','bold')
    set(    get(haxes(j),'YLabel')  ,'FontSize',14, 'FontName',fntName)   %,'FontWeight','bold')
    set(    get(haxes(j),'title')   ,'FontSize',14, 'FontName',fntName)    %,'FontWeight','bold')
    set(    haxes(j)                ,'FontSize',14, 'FontName',fntName)
    set(    haxes(j)                ,'LineWidth',1.0)

end

if nargin < 1
    % these are kind of arbitrary, but easy to change
    lbwh = [0.5871    2.7218    9.8198    7.9919];
end

%set(gcf,'PaperOrientation','Portrait')
set(gcf,'Units','Inches')         
set(gcf,'Position',lbwh)
set(gcf,'PaperPosition',lbwh)

% Weird bug in axes.m
% axes(htop);
% Try this:
set(gcf,'CurrentAxes',htop);


end

function publicationStandards_older(lbwh)
% PUBLICATION STANDARDS.M
% publicationStandards(lbwh)
% A comprehensive list of things to change 
% to make pub quality figs.
%
% Optionally input figure position (inches) 

% Copyright (C) 2000-2010 Brian Emery

% set(hfig,'DefaultAxesFontSize',14)
% set(hfig,'DefaultAxesLineWidth',1.5) 

fnt = 14;
fntName = 'TimesNewRoman';

% keep handle of current axes
htop = gca;

% Get all axes handles
haxes = findobj(get(gcf,'Children'),'type','axes'); 

% get rid of colorbars and legends
haxes = haxes(~strcmp('Colorbar',(get(haxes,'Tag'))));
haxes = haxes(~strcmp('legend',(get(haxes,'Tag'))));

% Loop over each axes and change settings
for j=1:length(haxes)
    
    % Removing this keeps all the layers in the right order
    % axes(haxes(j))
    
    set(    findobj('LineWidth',.5) ,'LineWidth',1.5)
    set(    findobj('FontSize',10)  ,'FontSize',fnt, 'FontName',fntName)
    set(    get(haxes(j),'XLabel')  ,'FontSize',14 , 'FontName',fntName)
    set(    get(haxes(j),'YLabel')  ,'FontSize',fnt, 'FontName',fntName)
    set(    get(haxes(j),'title')   ,'FontSize',fnt, 'FontName',fntName)
    set(    haxes(j)                ,'FontSize',fnt, 'FontName',fntName)
    set(    haxes(j)                ,'LineWidth',1.5)

end

if nargin < 1
    % these are kind of arbitrary, but easy to change
    lbwh = [0.5778    0.3444   12.2222    7.9556];
end

set(gcf,'PaperOrientation','Portrait')
set(gcf,'Units','Inches')         
set(gcf,'Position',lbwh)
set(gcf,'PaperPosition',lbwh)

axes(htop)


end

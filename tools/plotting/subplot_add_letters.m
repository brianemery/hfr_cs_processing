function subplot_add_letters(ax,x,y)
% SUBPLOT ADD LETTERS - add letters to subplots 
%  subplot_add_letters(ax, x, y)
%
% INPUT
% ax - array of axes handles
%
% OPTIONAL INPUTS
% % Set position relative to axes (these are default)
% x = 0.05;
% y = 0.9 ;
%
% EXAMPLE CUSTOMIZATION
% % reposition the letter 'C'
% h = findobj('String','C');
% pos = get(h,'Position');
% set(h,'Position',pos + [0 -25 0])


% Copyright(C) 2011 Brian Emery

% TO DO
% figure out how to pin this to the axes so that when they change the
% letters float with the axes (see annotations)

% abort if only one axis
if length(ax) == 1, return, end

% gca for use later
h1 = gca;

% get alphabet letters (lowercase)
letter = char(97:97+25);

if nargin < 2
    % Set position relative to axes
    x = 0.05;
    y = 0.9 ;
    
end

% get axes if not input
if nargin < 1
    ax = find_axes;
end


for i = 1:numel(ax)
    
    % Set current axes
    set(gcf,'CurrentAxes',ax(i))
    
    % Convert to data units
    a = axis(ax(i));
    
    xd = a(1) + (( a(2)-a(1) )*x);
    yd = a(3) + (( a(4)-a(3) )*y);
    
    
    % add text
    ht = text(xd,yd,upper(letter(i))); % ')']) %,'units','normalized') 
    
    set(ht,'FontName','Times New Roman') %,'FontSize',12) % <-- custom line for CF
    
end

set(gcf,'CurrentAxes',h1)

end

function ax = find_axes
% FIND AXES

% Get all axes handles
ax = findobj(get(gcf,'Children'),'type','axes'); 

% get rid of colorbars and legends
ax = ax(~strcmp('Colorbar',(get(ax,'Tag'))));
ax = ax(~strcmp('legend',(get(ax,'Tag'))));

ax = sort(ax);

end
function colorbar_label(str,h)
% COLORBAR LABEL.M
% colorbar_label(str,h)
%
% INPUTS
% str - string to use as label
% h   - optional colobar handle

% Copyright (C) 2019 Brian Emery 

% Notes about changeing colorbar labels:
% 
% h=colorbar;
%    set(get(h,'Ylabel'),'String','^oC','Rotation',0)

if nargin < 2
    h = findobj('Tag','Colorbar');
end

 % old way
 % set(get(y,'Ylabel'),'String','WindU Variance (cm^2s^-^2)')
 h.Label.String = str;
 

 
end
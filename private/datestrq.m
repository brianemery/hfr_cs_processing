function S = datestrq(D)
% convert a date number to a string with format yyyy_mm_dd_HHMMSS
%   datestrq(734474.62)
%       ans = 2010_12_02_145248
%  


% Copyright: zhang@zhiqiang.org, 2010
% Customized Nov 2010 Brian Emery ('yyyy_mm_dd_HHMM')
%   renamed datestrq.m for better recall

% TO DO
% take advantage of sprintf vectorization by converting cell to char?
% allow typical datestr inputs, then parse to create sprintf string

if iscell(D)
    S = cell(size(D));
    for i = 1:numel(D)
        S{i} = mdatestr(D{i});
    end
    return;
end

if ~ischar(D)
    if min(D) > 1000 && max(D) < 693960
        D = D + 693960;
    end
    
    D = datevecmx(D); 
    % yyyy_mm_dd_HHMMSS
    % S = sprintf('%04d_%02d_%02d_%02d%02d%', D(1), D(2), D(3), D(4), D(5));
    S = sprintf('%4.0f_%02.0f_%02.0f_%02.0f%02.0f%02.0f', D(1), D(2), D(3), D(4), D(5), D(6));
    
elseif numel(D) == 10 && D(8) == '-'
    S = D;
else
    S = formatdate(D);
end

end
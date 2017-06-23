function legh = legend_append(h,cel,opt)
% LEGEND APPEND - append to existing legend
% legh = legend_append(h,cel,opt)
%
% h   = new plot object handles
% cel = cell array of assoc labels
% opt = 1 to add to top of legend (defaults to 0)
%
% based on addTextLikeLegend

% Copyright (C) 2010 Brian Emery
% 4Feb2010  

if nargin<3, opt=0; end

if ischar(cel), cel = cellstr(cel); end

% get info if there is a pre-existing legend: (outm is cell)
[legh,objh,outh,outm] = legend;

% if no legend, make outm a cell
if isempty(outm), outm ={}; end

% create legend inputs
if opt
    hh=[h(:); outh];
    inArray={cel{:} outm{:}};
else
    hh=[outh; h(:)];
    inArray={outm{:} cel{:}};
end

legh=legend(hh,inArray);


end
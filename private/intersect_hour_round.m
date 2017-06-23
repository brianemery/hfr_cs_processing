function [c,ia,ib] = intersect_hour_round(A,B,hr)
% INTERSECT HOUR ROUND - intersect with rounding in hours
% [c,ia,ib] = intersect_hour_round(A,B,hr)
%
% Simplifies calls such as:
%
% [c,ia,ib] = intersect(round(A*24)/24,round(B*24)/24);
%
% Which round to the nearest hour, improving the use of 
% intersect.m
%
% hr input defaults to 24 if not provided. For nearest minute use hr =
% 1440, nearest 10 minutes hr = 144for example

% Copyright (C) 2011 Brian Emery

if nargin < 3, hr = 24; end

[c,ia,ib] = intersect(round(A.*hr)./hr,round(B.*hr)./hr);

end
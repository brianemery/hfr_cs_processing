function t = strjust(s,justify)
%STRJUST Justify character array.
%   T = STRJUST(S) or T = STRJUST(S,'right') returns a right justified 
%   version of the character array S.
%
%   T = STRJUST(S,'left') returns a left justified version of S.
%
%   T = STRJUST(S,'center') returns a center justified version of S.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 1997/11/21 23:47:38 $

if nargin<2, justify = 'right'; end

if isempty(s), t = s; return, end

% Find non-pad characters
ch = (s ~= ' ' & s ~= 0);
[r,c] = find(ch);
[m,n] = size(s);

% Determine offset
switch justify
case 'right'
    spa = sparse(c,r,c);
    offset = full(max(spa))';
    offset = n-offset;
case 'left'
    [dum,offset] = max(ch,[],2);
    offset = 1 - offset;
case 'center'
    spa = sparse(c,r,c);
    offset1 = full(max(spa))';
    [dum,offset2] = max(ch,[],2);
    offset = floor((n - offset1 - offset2 + 1)/2);
end

% Apply offset to justify character array
newc = c + offset(r);
t = repmat(' ',m,n);
t(r + (newc-1)*m) = s(r + (c-1)*m);

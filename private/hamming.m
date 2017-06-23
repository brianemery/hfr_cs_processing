function d = hamming(n)
% HAMMING.M	- returns the N-point Hamming window.
% d = hamming(n)
%
% Input		: n = number
%
% Output	: w = vector
%
% Usage		: w = hamming (n)
%
% Comments	: allows also the call:  hamming(xx), taking only format from signal xx
%
% See also	: HAMMING2

% modification of original MATLAB (3.5)  file
% HGFei, 1990  

n = n(:); 
% [hig, wid] = size(n);
% 
% if  wid > 1;  n = wid; end; 

d = (.54 - .46*cos(2*pi*(0:n-1)'/(n-1))).';

% if  hig > 1;
% ww = w;  
%  for jj = 2 : hig;
%  w(jj,:) =  ww;
%  end; 
% end; 

end


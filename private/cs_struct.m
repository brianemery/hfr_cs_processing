function CS = cs_struct(N)
% CS STRUCT - create empty CS structure
%
% INPUT
% N  - number of elements of CS

% Standard format has the matricies with
% dimensions nfft x # range cells


[CS.Header, ...
 CS.antenna1Self, ...
 CS.antenna2Self, ...
 CS.antenna3Self, ...
 CS.antenna12CrossSp, ...
 CS.antenna13CrossSp, ...
 CS.antenna23CrossSp, ...
 CS.spectraQualNum, ...
 CS.rdCSErr] = deal([]);

CS.FileName = '';

CS.ProcessingSteps = {''};
CS.Units = 'volts^2';

if nargin > 0 
   
    CS(1:N) = deal(CS);
    
end


end
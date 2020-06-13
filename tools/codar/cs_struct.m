function CS = cs_struct(N,M)
% CS STRUCT - create empty CS structure
%
% INPUT
% N  - number of elements of CS
% 
% NEW 2nd input, allows field names for arbitrary arrays with M elements
% 
% DATA SPECIFICATION
% Standard format has the matricies with
% dimensions nfft x # range cells
%

[CS.Header] = deal([]);


if nargin < 2

[CS.antenna1Self, ...
 CS.antenna2Self, ...
 CS.antenna3Self, ...
 CS.antenna12CrossSp, ...
 CS.antenna13CrossSp, ...
 CS.antenna23CrossSp, ] = deal([]);


else

  % make fields
  CS = make_fields(M,CS);
    
        
end


[CS.spectraQualNum, ...
 CS.rdCSErr] = deal([]);


CS.FileName = '';

CS.ProcessingSteps = {''};
CS.Units = 'volts^2';

if nargin > 0 
   
    CS(1:N) = deal(CS);
    
end

CS.freqs =[];
CS.Vrad = [];


end

function CS = make_fields(M,CS)


fn = cs_make_field_names(M);


for n = 1:length(fn)
        
        % now calc auto and cross products
        CS.(fn{n}) = [];
        
end


end


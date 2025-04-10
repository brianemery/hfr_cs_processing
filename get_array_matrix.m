function A = get_array_matrix(APM)
% GET ARRAY MATRIX - Array matrix from SeaSondes APM struct
%  A = get_array_matrix(APM)
% 
% A is the Array Manifold Matrix (c.f. Van Trees 2002, pg 31), aka the
% collection of vectors that incorporate all of the spatial characteristics
% of the array. 
%
% DATA FORMAT DEFINITION
% - The array matrix A is #elements x #bearings (THIS IS OUTPUT)
% - Similar data in the APM struct is the reverse, #bearings x #elements
% (historical reasons)

% Copyright (C) 2016 Brian Emery


% 
% if isfield(APM,'A') && ~isempty(APM.A)
%     A = APM.A.';
%     
% else
    
    % Get the array matrix from the APM struct
    A(1,:) = APM.A13R + 1i*APM.A13I;
    A(2,:) = APM.A23R + 1i*APM.A23I;
    A(3,:) = APM.A33R + 1i*APM.A33I;
    
% end
end
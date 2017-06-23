function A = get_array_matrix(APM)
% GET ARRAY MATRIX - Array matrix from APM struct for SeaSondes
%  A = get_array_matrix(APM)
% 
% A is the Array Manifold Matrix (c.f. Van Trees 2002, pg 31), aka the
% collection of vectors that incorporate all of the spatial characteristics
% of the array. 

% Copyright (C) 2016 Brian Emery

% Get the array matrix from the APM struct
A(1,:) = APM.A13R + 1i*APM.A13I; 
A(2,:) = APM.A23R + 1i*APM.A23I; 
A(3,:) = APM.A33R + 1i*APM.A33I; 

end
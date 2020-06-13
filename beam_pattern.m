function B = beam_pattern(A,ix)
% COMPUTE BEAM PATTERN - see Van Trees, 2002
% B = beam_pattern(A,ix)
%
% INPUT
% A  - Array Manifold (Matrix)
% ix - Bearing index of (column) to use in forming the beam pattern
% 
% Given one input (A), this assumes uniform weighing. 
% 
% Note that (abs(B)).^2 gives
% the power as a function of theta, and thus that this 
% can be used to find the half power beamwidth
%
% Used to validate make_ula_pattern.m 
%
% See also beamwidth.m 

% Copyright(C) 2017 Brian Emery
%
% May 2017

if strcmp(A,'--t'), test_case, return, end


[M,~] = size(A);

if nargin < 2
    % compute the beam pattern with uniform weighting
    wt = (1/M)*ones(M,1);  % eqn 2.90
    
else
    % otherwise define the 'conventional beam pattern' as on pg 53 which
    % essentially steers the direction and gives the beam pattern in that
    % direction. With this B is determined by equation 2.124
    wt = (1/M)*A(:,ix); 
    
end

% Eqn 2.69
B = wt'*A;




end

function test_case

% theta relative to array axis
th = 0:180;
M = 10;

[A,dm] = make_ula_pattern(th,M,'vt02');

% get the beam pattern
B = beam_pattern(A); 

% Plot half the bottom figure 2.18, pg 46, which plots the magnitude of B
plot(th,abs(B)), title('Compare wtih VT02 Fig 2.18, pg 46')

hold on

% Same as broadsize steering
B = beam_pattern(A, find(th==90) ); 

plot(th,abs(B),'o')


B = beam_pattern(A, find(th==0) );

keyboard

end

% first try
function p = compute_beamform(A)
% SEE TF09 pg 6 and try to make the beam pattern fig. This is close
%
% ... something not quite right, not sure what. Their figure is much closer
% to a 10 element array
%
% NOTATION NOTE
% abs(vector) in TF09 is the 2-norm (that is the sum(abs(vector))

% see pg 5 for frequency domain beam forming
for i = 1:size(A,2), 
    W = A(:,i)/norm(A(:,i));
    p(i) = norm( W'*A(:,81) ).^2; 
end

end 

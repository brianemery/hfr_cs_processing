function A = make_ura_pattern(phi,M,N,th)
% MAKE URA PATTERN - make uniform rectangular array manifold vector
% A = make_ura_pattern(th,M)
% 
% Based on Van Trees, 2002, Ch 4. The uniform rectancular array with
% elements on a grid with lambda/2 spacing (in both x and y), is known as
% the standard rectangular array (SRA). Typically d <= lambda/2 (see pg
% 239). 
%
% INPUTS
% N     - number of elements in x direction, 
% M     - number of elements in y direction 
% phi   - steering angles (in degrees relative to x (eg ccwE))*
%
% OUTPUTS
% A     - complex array matrix (#elements x #bearings) = N*M x length(th)
%
% SEE ALSO
% make_ula_pattern.m, make_ideal_pattern.m, make_uca_pattern.m

% Copyright (C) 2017 Brian Emery
%
% May 2, 2017

% TO DO
% finish codeing, test functionality
% need something to check against

% NOTES
% - Also, uniform weighting is assumed ? 
% - the theta input is currently for a single value of theta



% Optional test case
if strcmp(phi,'--t'), test_case, return, end

% Enable optional vertical pattern - for validation mostly
if nargin < 4
    th = 90; % degrees from vertical
end

% Define dx/lambda and dy/lambda, that is, normalized by lambda
dx = 0.5;
dy = 0.5; 




% See Van Trees, 2002, pg 235, 240, 249

% Make grid of M and N, indexed from zero, M columns
% % this seems to put the phase center at a corner
% [m,n] = meshgrid(0:M-1,0:N-1);
%
% Try to put the phase center in the middle
[m,n] = meshgrid(  (0:M-1) - mean(0:M-1)  ,   (0:N-1) - mean(0:N-1)   );

% vectorize these
m = m(:);
n = n(:);

% now expand to size of phi
m = repmat(m,1,length(phi));
n = repmat(n,1,length(phi));

% make matrix out of phi
phi = repmat(phi(:)',size(m,1),1);


% eqn 4.2, 4.3 with dx defined in wavelengths
psi_x = 2.*pi.*dx.*sind(th).*cosd(phi);
psi_y = 2.*pi.*dy.*sind(th).*sind(phi);


% Compute the matrix of array manifold vectors (eqn 4.50, 4.53)
A = exp( 1i.* ( (n .* psi_x)  + (m .* psi_y) ) );



end

function test_case
% Need a good check from somewhere that I can reproduce ...

% Consistency Check:
% 4 element rectangular and UCA should be the same 
% 8 element UCA can be two URAs ofset 45 deg ...
% M or N = 1 URA is basically a ULA 
%

% vs ULA 
N = 4; M=1; %M = 4; N = 1;
th = -180:180; 

[L,~] = make_ula_pattern(th,N,'vt02');

R = make_ura_pattern(th,M,N);

ix = find(th ==0);

 pl = beam_pattern(L,ix);
 pr = beam_pattern(R,ix);
% pl = compute_beamform(L,ix);
% pr = compute_beamform(R,ix);

plot(th,pl,'-o'), hold on
plot(th,pr,'*')

title('ULA and a 1-row URA')

keyboard


% vs UCA
% the UCA radius in this case is sqrt(2)*d/2 where d is the URA lattice
% spaceing (0.5 here) 
%
% STATUS: not exactly correct b/c zero phase location is different ...

M = 2; N = 2;
th = -180:180; 

C = make_uca_pattern(th,M+N, sqrt(2)*(0.5)/2);

R = make_ura_pattern(th,M,N);

ix = find(th ==0);

pc = beam_pattern(C,ix);
% pc = compute_beamform(C,ix);
% pr = compute_beamform(R,ix);

% The UCA puts an element on the x axis, so rotate the URA to match
ix = find(th == 45);

pr = beam_pattern(R,ix);

figure
plot(th,pc,'-o'), hold on
plot(th-45,pr,'*')

title('4 element UCA (0^o( and URA with 45^o beam pattern)')



keyboard





% TEST PLOT 

phi = -180:179;

A = make_ura_pattern(phi,3,3,0);


% see pg 5 for frequency domain beam forming
for i = 1:size(A,2), bf(i) = abs(A(:,i)'*A(:,81)); end

plot(th, 20*log10(bf))

axis([-80 80 -10 30])

keyboard

% Use Phased Array Toolbox for comparison?
% 
% See help phase.URA
%
% URA object can do triangular lattice also

% Example 1:
%   Construct a 2x3 URA and plot its azimuth response. Assume the
%   operating frequency is 1 GHz and the wave propagation speed is 3e8
%   m/s.

ha = phased.URA([3 3]);
fc = 13e6; c = 3e8;
plotResponse(ha,fc,c,'RespCut','Az','Format','Polar');

% see what it looks like ...
viewArray(ha)

% array gain, very similar to array steering vector for uniform weights
hag = phased.ArrayGain('SensorArray',ha);

gn = step(hag,3e8,[phi;zeros(size(phi))]);

% Example 2:
%   Find the response of each element in the above array at the
%   boresight.

ang = [0;0];
resp = step(ha,fc,ang)







end

function p = compute_beamform(A,ix)
% SEE TF09 pg 6 and try to make the beam pattern fig. This is close
%
% ... something not quite right, not sure what. Their figure is much closer
% to a 10 element array
%
% NOTATION NOTE
% abs(vector) in TF09 is the 2-norm (that is the sum(abs(vector))
%
% ix is the index of theta_zero,

% see pg 5 for frequency domain beam forming
for i = 1:size(A,2), 
    W = A(:,i)/norm(A(:,i));
    p(i) = norm( W'*A(:,ix) ).^2; 
end

end 


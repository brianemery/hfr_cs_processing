function [A,dm] = make_ula_pattern(th,M,swch)
% MAKE ULA PATTERN - uniform linear array
% [A,dm] = make_ula_pattern(th,M,swch)
%
% Construct the matrix of antenna location vectors 
%
% Based on Tuncer and Friedlander, 2009 pg 33, and others (see swch
% definition below). The default has th defined as relative to the array
% normal.
%
% INPUTS
% M     - number of elements
% th    - steering angles in degrees (zero at array normal, unless
%         otherwise noted, also positive in counter-clockwise (math) 
%         convention 
% swch  - optional, with possible values:
%
%         'krim', th has zero along direction of array axis (Krim and
%         Viberg, 1996)
%
%         'xu', uses definition on pg 2564 of Xu and Buckley, 1992, has
%         zero degrees (th = 0) on the array normal
%
%         'moses' uses code from Randy Moses, th = 0 on the array normal
%         (http://www2.ece.ohio-state.edu/~randy/SAtext/
%
%         'vt02' has theta = 0 along array axis, and array phase center in
%         middle of array
%
% OUTPUTS
% A     - complex array matrix (#elements x #bearings)
% dm    - element distances from array center in wavelengths (ie d/lambda)
%
%
% Note that arrays of this sort have a null to null beamwidth (BW_nn in Van
% Trees, 2002) of BW = 4/M (radians)
%
% The output is the array response to a unit amplitude signal, and is known
% as the array manifold. A output here is equivalent to the steering matrix
% given by Krim and Viberg, 1996.
%
% SEE ALSO
% make_phased_array_pattern.m which is basically incomplete but uses the
% APM struct format ...

% Copyright (C) 2017 Brian Emery



if strcmp(th,'--t'), test_case, return, end

if nargin < 3, swch = ''; end

switch swch
    
    case 'krim'
        
        % KRIM AND VIBERG 1996 METHOD
        % equation 12. Similar to Stoica and Nehorai 1989
        % k = 2pi/lambda, d = lambda/2 = interelement spacing, thus kd = pi
        %
        % First element at the array origin, theta relative to array axis?
        
        % get number of bearings
        n = length(th);
        
        % expand this into matricies for use in A definition
        L = repmat(0:M-1,n,1);
        th = repmat(th.',1,M);
        
        
        % This is sideways, each row is a bearing, cols are antennes,
        % fixed by transpose
        A =  ( exp( -1i.*pi.*cosd(th).*L ) ).';
        
        % half wavelength spacing is assumed, dm in wavelengths
        dm = (0:M-1)/2;
        
    case 'xu'
        
        % See Xu and Buckley, 1992, just after eqn 24
                
        % this is the matrix if integers, same size as output 
        n = repmat((0:M-1)',1,length(th));

        % expand bearings to size of output  
        th = repmat(th(:)',M,1);
        
         
        % a(th) = [1 exp(-j*pi*sin(th)) ... exp(-j*9*pi*sin(th) ]';
        % A = [ones(1,size(th,2)); exp(-1i.*pi*n
        A =  exp( -1i.*pi.*n.*sind(th) );
        
        dm = [];
        
    case 'moses'
        
        % from uladata.m, with d = 0.5.
        A = exp(-2*pi*1i*(0.5)*(0:M-1).'*sin((th(:).')*pi/180));

        dm = [];
        
    case 'vt02'
        
        % Van Trees, 2002, defines the ULA with the phase center in the
        % middle, and theta relative to the array axis. See pg 40 (eqn
        % 2.73), and eqn 2.72, with their n, N, renamed m, M here.
        dm = (0:(M-1))' - (M-1)/2;
        
        % the indexing (n - (N-1)/2), expanded out ...
        m = repmat( dm , 1, length(th) );
        
        % expand bearings to size of output
        th = repmat(th(:)',M,1);
        
        % fix spacing to lamba/2, (see eqn 2.59)
        % psi = 2*pi/lambda * cos(th) *d  % if d ~= lambda/2
        psi = pi .* cosd(th);
        
        A = exp( 1i .* m.* psi);  
                        
    otherwise
        
        % TF09 METHOD
        %
        % DEFINE ARRAY DIMENSIONS
        % Compare with TF09 equation 1.27, gets the same answer (see below)
        % relative to lambda ie d/l, usually l/2 appart, making the total length
        % to be (M-1)*(l/2), column vector
        m = (M-1)/2;
        
        % element distances from array center in wavelengths
        dm = (1/2)*(-m:m)';
        
        % expand to size of A output
        d = repmat(dm,1,length(th));
        
        % expand bearings to size of output
        th = repmat(th(:)',M,1);
        
        % compute array matrix
        A = exp( 1i.*2.*pi.*d.*sind(th) );
        
        
end

return

% NOTES




% TF09 equation 1.27
M = 3;
m = 1:M; % element index

% Define the array aperture (total length in wavelengths)
L = (M-1)/2;

% compute element distances from array center
dm = (m-1)*(L/(M-1)) - (L/2)








end

function test_case
% TEST PLOT

% SEE ULADATA.M from Spectral Analysis of Signals by Stoica and Moses
%
% Y=uladata(theta,P,N,sig2,m,d);
%
%    theta  <- arrival angles of the m sources in degrees
%    P      <- The covariance matrix of the source signals
%    N      <- number of snapshots to generate
%    sig2   <- noise variance
%    m      <- number of sensors
%    d      <- sensor spacing in wavelengths 
%    Y      -> m x N data matrix Y = [y(1),...,y(N)]
%
% Copyright 1996 by R. Moses
% 
% % definitions
% th = -80:80;
% d = 0.5;
% M = 16;

% COMPARE WITH VAN TREES 2002
%
% Figure 2.18 

% theta relative to array axis
th = 0:180;
M = 10;

[A,dm] = make_ula_pattern(th,M,'vt02');

% get the beam pattern
B = beam_pattern(A); 

% Plot half the bottom figure 2.18, pg 46, which plots the magnitude of B
plot(th,abs(B)), title('Compare wtih VT02 Fig 2.18, pg 46')


% plot empirical and analytical

plot(th,B), hold on

t = th*pi/180;
Ba = (1/M) .* sin( (M/2) .* pi .* cos(t) )./sin( (1/2) .* pi .* cos(t) );

plot(th, Ba,'o')










keyboard

% Compare beam patterns
p1 = compute_beamform(A1);

p = compute_beamform(A);


plot(th, 10*log10(p), '-*'), hold on
plot(th, 10*log10(p1), '-o'), hold on

axis([-80 80 -60 30])

keyboard

% Here is a polar plot that uses windowing some how:
% http://ifmaxp1.ifm.uni-hamburg.de/WERA_Guide/Hamming-00.gif
% http://ifmaxp1.ifm.uni-hamburg.de/WERA_Guide/Hamming-45.gif


return


% Compare with MATLAB PHASED ARRAY TOOLBOX
%
% Status: toolbox is a PITA and not helpful. Cant tell what is going on ...

th = -180:180;
M = 8;

[A,~] = make_ula_pattern(th,M);

% see pg 5 for frequency domain beam forming, this 
for ii = 1:size(A,2), pwr(ii) = norm(A(:,ii)'*A(:,181)).^2; end

plot(th,20*log10(pwr/max(pwr)))
axis([-200 200 -80 0])

% MATLAB PHASED ARRAY TOOLBOX
% ... is extremely annoying. "I know, let's just hide everyting so the dumb
% user doesn't have to understand what is going on" -mathworks
%
fc = 13e6; % freq in MHz

hula = phased.ULA('NumElements',M,'ElementSpacing',0.5);

hsv = phased.SteeringVector('SensorArray',hula);

SV = step(hsv,fc,[th; zeros(size(th))]);




end




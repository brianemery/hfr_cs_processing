function [idx,Vsml] = sml_ap(A,R,n)
% SML VIA ALTERNATING PROJECTION -  Stocastic ML DOA with MLE-AP algorithm
% idx = sml_ap(A,C,q)
%
% Stocastic Maximum Likelihood per Stoica and Nehorai 1990 where it is 
% Know as the Unconditional ML, and Krim and Viberg 1996. 
%
% INPUTS
% A    - Array Manifold (aka array matrix at all theta (complex))
% C    - Covariance matrix*
% q    - max number of signal sources**
%
%  * Must be positive definite (invertible), usually true if 
%    K (snapshots) >= N (antennas)
% ** If a radar has M elements and theta has d bearings, A is Mxd
%
%
% OUTPUTS
% idx    - Cell array with indecies of the n signal source solutions
%          eg idx{1} has the index of the single source soln, idx{2} is
%          dual, etc
% Vsml   - The log of the determinant of the structured estimate of R, 
%          aka R^, output for comparison with R in the detection problem
%          (See glrt.m)
% 
%
% REFERENCE
% Krim and Viberg, 1996, "Two Decades of Array Signal Processing Research: 
%   the Parametric Approach", IEEE
% Van Trees, 2002, "Optimum Array Processing, Part IV", section 8.5.1
%
%
% SEE ALSO
% music.m, mle.m, mle_ap.m, glrt.m

% Copyright (C) 2017 Brian Emery
%
% Version 23-Jun-2017 15:59:26



% TO DO
% - test for accuracy/validity (simulation with reproduce_zw88.m)
% - output cell array of indecies like music.m?
% - does this favor the ends of the APM with real data?
% - check for 3 bearing solutions when > 3 antennas? 
%
% clean up arg_max code 

% DONE
% - functional test 


% check for test case
if strcmp('--t',A), test_case;  end



% INITIALIZATION - use MLE-AP for first pass 
% using eqn 17 and 18, etc (ZW88)
thi = [];
for i = 1:n
    
    % get indecies of the theta_i
    ix = arg_max(A,R,thi);
    
    % add them in to make the augmented matrix
    thi = [thi ix];
end



% MAIN LOOP
d = ones(size(n)); 
k = 1;
thi_k = thi; % thi is current iteration, thi_k is the next

while any(d > eps) 
    
    % at each iteration, loop over each possible signal, then compare
    % result with the previous until they stay the same
    
    for i = 1:n
        
        % index of thi to use in prior projection
        x =  setdiff(1:n,i);
                
        % loop over th, for  this iteration (note this n is not an index)
        [Vsml,thi_k(i)] = arg_min(A,R,thi(x),n);
        
                
        % disp(['thi   = ' num2str(thi)])
        % disp(['thi_k = ' num2str(thi_k)])
        
        % check for convergence for each individually   
        d = sqrt(( thi_k - thi ).^2);
        
        % update/ inti thi_k
        thi = thi_k;
        
    end
    
    
    k = k+1;  % disp(num2str(k))
    
    % Put a break in here  (Viberg et al 1991 uses max 30 iterations)
    if k >= 100, % see experiment_mle_ap_iterations
        disp('SML: 100 iterations reached, terminated')
        break
    end

    
    
end

idx = thi(:);


end

function [Vsml,thi] = arg_min(A,R,thx,M)
% ARG MIN - the function to minimize
%
% INPUTS
% A   - the manifold matrix (L rows = # of antennas, n bearings)
% R   - data covariance matrix
% thx - the index of the 'prior' solution to leave out of the search
% M   - the signal subspace dimension (aka number of signals)
% 
% Krim and Viberg 1996, equations 59, 60, 61
% Ottersten, Viberg, Stoica and Nehorai, 1993 eqns 4.37, 4.38 


% Get array elements and number of bearings to loop over 
[L,n] = size(A);

I = eye(L);


% inti storage
Vsml = NaN(n,1); 

for i = setdiff(1:n,thx)
    
    % get pseudo inverse of the part of A that we want
    At = pinv(A(:,[thx i]));
    
    % Projection matrix (function of theta) PI_A perp
    PI = I - (A(:,[thx i])*At);  % Pv_perp for VT02
   
    % sigma sml squared as a function of theta
    % Equation 58, theirs is missing an R
    sigma = trace(PI*R)/(L-M);
    
    % eqn 60 (This is S^ (S hat) in Ottersten et al 1993)
    P = At*(R-(sigma*I))*At'; 
    
    % eqn 61
    % compute function of which to find minimum
    % Stoica and Nehorai have this as the
    % ln of the determinant ... need to check (eqn 2.20)
    % "log(X) is the natural logarithm of the elements of X".
     Vsml(i) = log( det( A(:,[thx i])*P*A(:,[thx i])' +  (sigma*I) ) ); 
    
end

% % uncomment for testing to see the incrementing DOA function
% plot(real(Vsml)), hold on, pause

% essential to take real part first
[~,thi] = min(real(Vsml));

% output just the solution for the detection problem
Vsml = real(Vsml(thi));
 

end


% MLE-AP code for initialization 
function thi = arg_max(A,R,ix)
% INIT MLE -  (eqn 17)
%
% loop over the possible values of theta_i (just the index)

% if there is an empty 3rd input, use this to get the augmented matrix
% - do I need to leave out an i? ... not in this part ... just init 
if nargin > 2
    B = A(:,ix);
else
    B = []; ix = [];
end


% Get number of bearings to loop over for the maximization
n = size(A,2);

% init storage
tr = NaN(n,1);

% loop over index of possible, excluding 
for i = setdiff(1:n,ix)
    
    P = proj([B A(:,i)]);
    
    tr(i) = trace(P*R);
    
end

[~,thi] = max(real(tr)); % need to verify?

end

function P = proj(A)
% COMPUTE PROJECTION MATRIX -
%
% Note: ' is the Hermitian transpose, which is intended

% P = A*(inv(A'*A)*A'); % can I remove the inverse?   % < --- SLOW
P = A*( (A'*A)\A' ); % inverse removed

end


% dev testing
function thi = arg_min_dev(A,R,thx,M)
% ARG MIN - the function to minimize - development version 
%
% ... same as above but loops over all the bearings to see the DOA function
%
% INPUTS
% A   - the manifold matrix (L rows = # of antennas, n bearings)
% R   - data covariance matrix
% thx - the index of the 'prior' solution to leave out of the search
% M   - the signal subspace dimension (aka number of signals)
% 
% Krim and Viberg 1996, equations 59, 60, 61


% Get array elements and number of bearings to loop over 
[L,n] = size(A);

I = eye(L);



% inti storage
tr = NaN(n,1);  trv = tr;

for i = 1:n % n here is total number of bearings to loop over
    
    % get pseudo inverse of the part of A that we want
    At = pinv(A(:,[thx i]));
    
    % Projection matrix (function of theta) PI_A perp
    PI = I - (A(:,[thx i])*At);  % Pv_perp for VT02

    
    % sigma sml squared as a function of theta
    % Equation 58, theirs is missing an R
    sigma = trace(PI*R)/(L-M);
    
    % eqn 60
    P = At*(R-(sigma*I))*At'; 
    
    % eqn 61
    % compute function of which to find minimum
    % Stoica and Nehorai have this as the
    % ln of the determinant ... need to check (eqn 2.20)
    % "log(X) is the natural logarithm of the elements of X".
    tr(i) = log( det( A(:,[thx i])*P*A(:,[thx i])' +  (sigma*I) ) ); 
%     
%     
    % VT02 has this instead .... which is identical but is missing the
    % the log . 
    %
    Pv =  A(:,[thx i])*pinv(A(:,[thx i]));
    
    % define Pv perp
    Pvp = I - Pv;
    
    % eqn 8.315 ... actually
    trv(i) = det( Pv*R*Pv +  ( trace(Pvp*R)*Pvp/(L-M) )  );
    

    
end

keyboard

% essential to take real part first
[~,thi] = min(real(tr));

 plot(real(tr)), hold on


 
end


function test_case
% TEST CASES
%
% NEW MLE-AP.M ...
% Elapsed time is 0.622195 seconds.
% Old mle_ap.m ...
% Elapsed time is 1.641843 seconds.
% mle_ap.m ...
% Elapsed time is 1.483485 seconds.
% mle.m ...
% Elapsed time is 228.815552 seconds.
% mle vs mle_ap ... get same result!
%
% SML-AP NOTE
% if M=N (elements and signals same) this blows up because of the sigma
% calculation which divides by #elements - #signals



% BASIC VALIDATION TEST
% high SNR, high K, ... works

th = -45:0.1:45;
M = 8;          % # elements in array
K = 100;        % # snapshots

% Make the pattern (xu seems to work best?)
[A,~] = make_ula_pattern(th,M,'xu'); %,'moses'); or 'xu', all relative to perp

[~,ix,~] = intersect(th,[3 14 25]);


% SNR 30
[C, asnr] = radar_simulation_basic( A(:,ix) ,K,40);


iy = sml_ap(A,C,3);


a = axis;
plot(ix([1 1]),a([3 4]))
plot(ix([2 2]),a([3 4]))

keyboard

dx = mle_ap(A,C,2);




% VALIDITY TEST
% From VT02, pg 991

reproduce_vt02_sml_ap



% DEVELOPMENT TEST
% look at the DOA function for simplest case

th = -45:1:45;
M = 3;          % # elements in array
K = 10;        % # snapshots

% Make the pattern (xu seems to work best?)
[A,~] = make_ula_pattern(th,M,'xu'); %,'moses'); or 'xu', all relative to perp

[~,ix,~] = intersect(th,3);


% SNR 30
[C, asnr] = radar_simulation_basic( A(:,ix) ,K,30);


thi = arg_min_dev(A,C,[],1);






end




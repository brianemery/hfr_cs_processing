function [idx,Vwsf,s] = wsf(A,R,n)
% WEIGHTED SUBSPACE FITTING -  Estimate DOA with WSF algorithm
% [idx,Vwsf,sigma2] = wsf(A,C,q)
%
% See Van Trees 2002, pg 1011, which suggests this is equivalent (for our
% purposes) to the Method of Direction Estimation (MODE) given by Stoica
% and Sharman, (1990).
%
% INPUTS
% A    - Array Manifold (aka array matrix at all theta (complex))
% C    - Covariance matrix
% n    - max number of signal sources (defaults to M-1)*
%
% *If a radar has M elements and theta has d bearings, A is Mxd
%
%
% OUTPUTS
% idx    - double array with indecies of the n signal source solutions
% Vwsf   - value of the the WSF cost function at the DOA estimate (theta)
% sigma2 - the estimated noise variance
%
% REFERENCE
% Krim and Viberg, 1996, "Two Decades of Array Signal Processing Research: 
%   the Parametric Approach", IEEE
% Viberg, Ottersten and Kailath, 1991, "Detection and Estimation in Sensor 
%   Arrays Using Weighted Subspace Fitting", IEEE  
%
% SEE ALSO
% music.m, mle.m, mle_ap.m, sml_ap.m

% Copyright (C) 2017 Brian Emery
%
% Version 13-Jun-2017 


% TO DO
% - test for accuracy/validity (simulation with reproduce_zw88.m)
% - output cell array of indecies like music.m?
% - does this favor the ends of the APM with real data?
% - check for 3 bearing solutions when > 3 antennas
%
% clean up arg_max code 

% DONE
% - eig "all eigenvectors normalized to unit norm"
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



% INIT WSF 

% Do eigen decomposition of covariance matrix input (R2013 compatible)
[U,D] = eig(R); D = diag(D); %,'vector'); 

% sort largest to smallest 
[~,ix] = sort(D,'descend'); 
D = D(ix);
U = U(:,ix); 

% Noise and Signal sub spaces for n brg soln
Us = U(:,1:n);
% Un = U(:,n+1:end);

% Matrix of Signal eigenvalues
Ds = diag(D(1:n)); 

% estimate noise variance from noise eigen values (aka sigma squared)
s = mean(D(n+1:end));


% Compute the Weight matrix (optimal for large K)
% VT02 eqn 8.368
W = (( Ds - s.*eye(n) )^2) /(Ds); % /(Ds) is eqivalent to *inv(Ds)




% MAIN LOOP
d = ones(size(n)); 
k = 1;
thi_k = thi; % thi is current iteration, thi_k is the next

while any(d > eps) 
    
    
    for i = 1:n
        
        % index of thi to use in prior projection
        x =  setdiff(1:n,i);
                        
        % loop over th, for  this iteration
        [thi_k(i),Vwsf] = arg_min(A,Us,W,thi(x));
                
        % check for convergence for each individually   
        d = sqrt(( thi_k - thi ).^2);
        
        % update/ inti thi_k
        thi = thi_k;
        
        
    end
    
    k = k+1;  % disp(num2str(k))
    
    % Put a break in here  (Viberg et al 1991 uses max 30 iterations
    if k >= 100, % see experiment_mle_ap_iterations
        disp('WSF: 100 iterations reached, terminated')
        break
    end
    
    
end

idx = thi(:);




end

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


[~,thi] = max(real(tr));

end

function [thi,Vwsf] = arg_min(A,Us,W,thx)
%  WSF Function to minimize
% 
% INPUTS
% A   - array matrix
% Us  - signal eigenvectors
% W   - weight matrix
% thx - index (?) of bearing to hold fixed for this iteration


% Get number of bearings to loop over for the maximization
n = size(A,2);

I = eye(size(A,1));


% inti storage
tr = NaN(n,1);

for i = setdiff(1:n,thx)
    
    % Projection matrix (function of theta) PI_A perp
    PI = I - (A(:,[thx i])*pinv(A(:,[thx i])));
    
    % compute trace
    tr(i) = trace( PI*Us*W*Us' ); 
    
    
end

% note that they dont take the log but sml_ap does and this looks much
% better (more accurate?) when taking the log
[~,thi] = min(real(log(tr)));

% plot(real(log(tr))), hold on
% 
% plot(thi,real(log(tr(thi))),'*'), pause

Vwsf = real(tr(thi));



end

function P = proj(A)
% COMPUTE PROJECTION MATRIX -
%
% Note: ' is the Hermitian transpose, which is intended

% P = A*(inv(A'*A)*A'); % can I remove the inverse?   % < --- SLOW
P = A*( (A'*A)\A' ); % inverse removed

end



function test_case
% DEV TEST
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



% BASIC VALIDATION TEST
% high SNR, high K, ... works

th = -45:0.1:45;
M = 8;          % # elements in array
K = 100;        % # snapshots

% Make the pattern (xu seems to work best?)
[A,~] = make_ula_pattern(th,M,'xu'); %,'moses'); or 'xu', all relative to perp

[~,ix,~] = intersect(th,[3 14 33]);


% SNR 30
[C, asnr] = radar_simulation_basic( A(:,ix) ,K,40);


iy = wsf(A,C,2);


a = axis;
plot(ix([1 1]),a([3 4]))
plot(ix([2 2]),a([3 4]))

keyboard







%   205   330
[C,APM] = test_data_tonys;

A = get_array_matrix(APM);


% NEW TEST
tic
ix = wsf(A,C,2);
toc


keyboard

dx = mle_ap(A,C,2);

end





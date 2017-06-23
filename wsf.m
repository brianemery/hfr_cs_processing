function idx = wsf(A,R,n)
% WEIGHTED SUBSPACE FITTING -  Estimate DOA with WSF algorithm
% idx = wsf(A,C,q)
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
% idx    - Cell array with indecies of the n signal source solutions
%          eg idx{1} has the index of the single source soln, idx{2} is
%          dual, etc
%
% REFERENCE
% Krim and Viberg, 1996, "Two Decades of Array Signal Processing Research: 
%   the Parametric Approach", IEEE
% Viberg, Ottersten and Kailath, 1991, "Detection and Estimation in Sensor 
%   Arrays Using Weighted Subspace Fitting", IEEE  
%
% SEE ALSO
% music.m, mle.m, mle_ap_ss.m, mle_ap.m

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

% Do eigen decomposition of covariance matrix input
[U,D] = eig(R,'vector'); 

% sort largest to smallest 
[~,ix] = sort(D,'descend');  % < --- CHECK
D = D(ix);
U = U(:,ix); 

% Noise and Signal sub spaces for n brg soln
Us = U(:,1:n);
% Un = U(:,n+1:end);

% Matrix of Signal eigenvalues
Ds = diag(D(1:n)); 

% estimate noise variance from noise eigen values % < --- CHECK
% does this need to be squared?
s = mean(D(n+1:end));


% Compute the Weight matrix
W = (( Ds - s.*eye(n) )^2) /(Ds); % /(Ds) is eqivalent to *inv(Ds)

% % what size I
% I = eye(m); %?









% MAIN LOOP
d = ones(size(n)); 
k = 1;
thi_k = thi; % thi is current iteration, thi_k is the next

while any(d > eps) 
    
    
    for i = 1:n
        
        % index of thi to use in prior projection
        x =  setdiff(1:n,i);
        
                
        % loop over th, for  this iteration
        thi_k(i) = arg_min(A,Us,W,thi(x));
        
                
        % disp(['thi   = ' num2str(thi)])
        % disp(['thi_k = ' num2str(thi_k)])
        
        % check for convergence for each individually   
        d = sqrt(( thi_k - thi ).^2);
        
        % update/ inti thi_k
        thi = thi_k;
        
        k = k+1;  % disp(num2str(k))
        
    end
    
    % Put a break in here?
    if k == 100,
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

[~,thi] = max(tr);

end

function thi = arg_min(A,Us,W,thx)
% 
% Simple version of MLE-AP ... for use with WSF.M


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

[~,thi] = min(tr);

% plot(real(tr)), hold on

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


% delete?



function idx = mle_ap(A,R,q)
% MLE VIA ALTERNATING PROJECTION -  Estimate DOA with MLE-AP algorithm
% idx = mle_ap(A,C,q)
%
% INPUTS
% A    - Array Manifold (aka array matrix at all theta (complex))
% C    - Covariance matrix
% q    - max number of signal sources (defaults to M-1)*
%
% *If a radar has M elements and theta has d bearings, A is Mxd
%
%
% OUTPUTS
% idx    - Cell array with indecies of the n signal source solutions
%          eg idx{1} has the index of the single source soln, idx{2} is
%          dual, etc
%
%
% SEE ALSO
% music.m, mle.m, mle_ap_ss.m

% Copyright (C) 2017 Brian Emery
%
% Version 22-May-2017 12:21:42


% TO DO
% - output cell array of indecies like music.m?
% - test, test and test 
% - does this favor the ends of the APM with real data?

% check for test case
if strcmp('--t',A), test_case;  end


% INITIALIZATION
% using eqn 17 and 18, etc
thi = [];
for i = 1:q
    
    % get indecies of the theta_i
    ix = arg_max(A,R,thi);
    
    % add them in to make the augmented matrix
    thi = [thi ix];
end


% MAIN LOOP
d = ones(size(q)); 
k = 1;
thi_k = thi; % thi is current iteration, thi_k is the next

while any(d > eps) 
    
    
    for i = 1:q
        
        % index of thi to use in prior projection
        x =  setdiff(1:q,i);
        
        % compute the projection for the non-i soureces
        Pb = proj(A(:, thi(x) ));
        
        
        % loop over th, for  this iteration
        thi_k(i) = arg_max_b(A,R,Pb,thi(x));
        
                
        % disp(['thi   = ' num2str(thi)])
        % disp(['thi_k = ' num2str(thi_k)])
        
        % check for convergence
        d = sqrt(( thi_k - thi ).^2);
        
        % update/ inti thi_k
        thi = thi_k;
        
        k = k+1;  % disp(num2str(k))
        
    end
    
    % Put a break in here?
    if k == 100,
        disp('MLE: 100 iterations reached, terminated')
        break
    end

    
end

idx = thi(:);

end

function thi = arg_max(A,R,ix)
% INIT MLE -  (eqn 17)
%
% loop over the possible values of theta_i (just the index)

% if there is a 3rd input, use this to get the augmented matrix
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

function P = proj(A)
% COMPUTE PROJECTION MATRIX -
%
% Note: ' is the Hermitian transpose, which is intended

% P = A*(inv(A'*A)*A'); % can I remove the inverse?   % < --- SLOW
P = A*( (A'*A)\A' ); % inverse removed

end

function b = proj_update(C,Pb)
% compute the b unit vector given/used in eqn 22
%
% given the update (a(theta)) and the prior projection (P_A for example) and
%
% USAGE
% b = proj_update(A(:,i),P)

I = eye(size(Pb));

% 'a(theta) sub A(THETA)'
Cb = (I - Pb)*C; 

b = Cb/norm(Cb); % eqn 22.b



end

function thi = arg_max_b(A,R,Pb,thx)
% b projection case (eqn 22.a)
%
% Form the augmentend matrix for multi-bearing cases


% Get number of bearings to loop over for the maximization
n = size(A,2);

% inti storage
tr = NaN(n,1);

for i = setdiff(1:n,thx)
    
    b = proj_update(A(:,i),Pb);
    
    tr(i) = b'*R*b;
    
end

[~,thi] = max(tr);

% plot(real(tr)), hold on

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
disp('NEW MLE-AP.M ...')
tic
ix = mle_ap(A,C,2);
toc


% DEV TEST CONTINUED

disp('Old mle_ap.m ...')
tic
% OLD MLE AP FOR COMPARISON
[oix,OAP] = mle_ap_ss(APM,C);
toc

disp('mle_ap.m ...')
tic
% MLE AP FOR COMPARISON
[absix,AP] = mle_ap_older(A,C,APM.BEAR);
toc

% keyboard
% 
% disp('mle.m ...')
% tic
% % MLE FULL SEARCH
% [ix,TR] = mle(APM,C);
% toc
% 





if isequal(ix,sort(absix(2:3)))
    disp('new mle_ap vs old mle_ap ... get same result!')
else
    disp('new mle_ap vs old mle_ap ... INCONSISTENT!')
end

keyboard

% PLOT THE 2 BRG SEARCH SURFACE
% 
% older code - keep for reference

[X,Y] = meshgrid(APM.BEAR,APM.BEAR);

h = surf(X,Y,real(TR)); set(h,'EdgeColor','Interp')

[mx,r,c] = max_2d(TR);

view(2)

z = get(h,'ZData');

set(h,'ZData',z-10)


% PLOT MLE AP SEARCH LINES

hold on

h = line(AP.x,AP.y); %real(AP.maxtr))


[DOA,singleIdx,dualIdx] = music(APM,C);

keyboard


end






% PREVIOUS VERSION
% - minor errors including: looking at the trace rather than bearings for
% convergence
% - two bearings search only
function  [absix,TR] = mle_ap_older(A,C,th) %,n) % APM,C)
% MLE VIA ALTERNATING PROJECTION -  Estimate DOA with MLE-AP algorithm
% [idx,TR] = mle_ap(A,C,th,n)
%
% INPUTS
% A    - Array Manifold (aka array matrix at all theta (complex))
% C    - Covariance matrix
% th   - bearings associated with A
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
%
% SEE ALSO
% music.m, mle.m, mle_ap_ss.m

% Copyright (C) 2017 Brian Emery
%
% Version 25-Apr-2017 12:21:42



% TO DO
% - Generalize for n sources (Up to n = M?)
% - output cell array of indecies like music.m?
% - Combine the subfunctions into one, pass indecies only?

% DONE
% Generalize for use with antenna array matrix




% SINGLE BEARING SEARCH
% Find the single bearing solution
% will probably need the value of the trace to compare with dual brg
% solution
[tr,ix] = find_mle(A,C);



% TWO BEARING SEARCH
% now, combine this (one brg) soln with the other elements of the apm,
% compute and store the trace

% First compute the location of the 2nd source, assuming the first is at ix
% from above

[absix,maxtr] = deal(NaN(15,1));

% Initialize the (absolute) bearing index, and the values of the trace
absix(1) = ix;
maxtr(1) = tr;

% INIT THE ITERATIVE ALGORITHM
m = 2;
dtr = 1;

% for m = 2:10
while dtr > eps
    
    % find the abs index of the APM to fix
    [tr,ix] = find_mle_dual(A,C,ix);
    
    
    maxtr(m) = tr;
    absix(m) = ix;
    
    % check for convergence
    dtr = abs( real(maxtr(m)) -  real(maxtr(m-1)) );
    
    
    % Put a break in here?
    if m == 25,
        % disp('MLE: 25 iterations reached, terminated early')
        break
    end
    
    % increment m
    m = m + 1;
    
    
end



% CREATE OUTPUTS
% output the indecies of the single and dual solutions


% strip off nan's
absix = absix(~isnan(absix));
maxtr = maxtr(~isnan(maxtr));


% Create the data for plotting the 2-D search line
TR.absix = absix;
TR.maxtr = maxtr;


% create the ordered pairs of 2 bearing solutions
TR.x = th([1; absix(2:2:end)]);
TR.y = th(absix(1:2:end));


if length(absix) > 3
    
    % grab the three (1 single and 2 dual) bearing solutions
    absix = absix([1 end-1 end]);
    
elseif length(absix) < 3
    % if we get here, somethign broken
    keyboard
end


return

keyboard

% PLOT TO LOOK FOR CONVERGENCE
plot(1:numel(maxtr),real(maxtr),'-o')


% keyboard
%
% [tr3,ix3,b3]= find_mle_dual(APM,C,ix2);
%
keyboard
%

% SINGLE VS DUAL (VS 3?) HYPOTHESIS TESTING



end


function [mxtr,ix] = find_mle(A,C)
% FIND MLE - single bearing case

% Get number of bearings to loop over
nt = size(A,2);

tr = NaN(nt,1);

for n = 1:nt
    
    % PRE CONSTRUCT THE APM MATRIX
    % A = [ APM.A13R(n)+ 1i.*APM.A13I(n)  APM.A23R(n) + 1i.*APM.A23I(n)  1+1i.*0 ].';
    
    % Note: ' is the Hermitian transpose, which is intended
    
    % P = A*(inv(A'*A)*A'); % can I remove the inverse?   % < --- SLOW
    %     P = A*( (A'*A)\A' ); % inverse removed
    
    %    P = A(:,n) * (inv( A(:,n)' * A(:,n) )) * A(:,n)'; % can I remove the inverse?
    P = A(:,n) * ( ( A(:,n)' * A(:,n) ) \ A(:,n)');
    
    tr(n) = trace(P*C);
    
    
end

% Ok ... what have I got?
[mxtr,ix] = max(tr);

end

function [mxtr,ix] = find_mle_dual(A,C,ii)
% FIND MLE - loop over APM bearings, finding the best fit to the signal

nt = size(A,2);

% create outputs
tr = NaN(nt,1);

% just do this with indexing, loop over index excluding ii


for n = setdiff(1:nt,ii)
    
    
    % Here A is A(theta), and P is a function of theta (thus loop over theta?)
    %
    % Note: ' is the Hermitian transpose, which is intended
    % This is eqn 11.b in Ziskind and Wax 1988
    % P = A*(inv(A'*A)*A'); % can I remove the inverse?   % < --- SLOW
    %     P = A*( (A'*A)\A' ); % inverse removed
    P = A(:,[ii n]) * ( ( A(:,[ii n])'*A(:,[ii n]) )\A(:,[ii n])' ); % inverse removed
    
    tr(n) = trace(P*C); % eqn 13.a  % < --- SLOW, maybe a faster way is tto just compute diag elements by hand?
    
    
end


% Ok ... what have I got?
[mxtr,ix] = max(tr);


% % % TEST PLOT
% % figure
% % plot(B.BEAR, real(tr)), hold on,
% % plot(B.BEAR(iy),real(tr(iy)),'r*')
%
%
% % ... so my two solutions so far are ix, and ii
% % need to translate ix into the APM index
%
% keyboard
%
% % TRANSLATE IX TO ABS INDEX
% ix = find(APM.BEAR == B.BEAR(iy));
%
%
% TO DO
% look at A, P at the two bearings with really high tr and make sure the
% math is working right ...

end

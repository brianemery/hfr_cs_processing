function  [absix,TR] = mle_ap_ss(APM,C)
% MLE VIA ALTERNATING PROJECTION - single and dual DOA estimates via MLE AP
% absix = mle_ap(APM,C)
%
% -- method for direction of arrival in SeaSonde processing --
% Was mle_ap.m, which is now being generalized. This here for bw
% compatibility
%
%
% INPUTS
% APM  - Standard APM struct
% C    - data covariance matrix (from cross spectra struct)
%
% OUTPUT
% absix - absolute index of the APM [single dual dual] at the estimated
%         direction of arrival
%


% TO DO
% Generalize for use with antenna array matrix
%
% Output whatever is needed to plot the 2 bearing map of the search 
%
% Test for single/dual? 
% 
% probably efficiencies to be had pre-computing all the P matricies first
% and storing them?
%
% Other improvements based on refs that cite this method?
% More extensive study?
%
% MLE/MUSIC outline:
% load/get CS
% get first order region (whole CS or range cell?)
% run doa (loop over range cells?)
% run single dual tests
% make the radial file 

% tic

% check for test case
if strcmp('--t',APM), test_case;  end

% SINGLE BEARING SEARCH
% Find the single bearing solution
% will probably need the value of the trace to compare with dual brg
% solution
[tr,ix] = find_mle(APM,C);



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
    [tr,ix] = find_mle_dual(APM,C,ix);
    
    
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

% toc





% CREATE OUTPUTS
% output the indecies of the single and dual solutions


% strip off nan's
absix = absix(~isnan(absix));
maxtr = maxtr(~isnan(maxtr));


% Create the data for plotting the 2-D search line
TR.absix = absix;
TR.maxtr = maxtr;


% create the ordered pairs of 2 bearing solutions
TR.x = APM.BEAR([1; absix(2:2:end)]);
TR.y = APM.BEAR([absix(1:2:end)]);


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


function [mxtr,ix] = find_mle(APM,C)
% FIND MLE - single bearing case

% Get number of bearings to loop over
nt = length(APM.BEAR);

tr = NaN(nt,1);

for n = 1:nt
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    % make it a column vector (M element x d signals), M = 3
    % A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
    
    % PRE CONSTRUCT THE APM MATRIX
    % If nargin == 2 we're looking for a 2 brg solution and need to modify the
    % APM search    
    A = [ APM.A13R(n)+ 1i.*APM.A13I(n)  APM.A23R(n) + 1i.*APM.A23I(n)  1+1i.*0 ].';

    
    
    % max(trace(P*C))
    % C is the covariance matrix (R in the //Introduction to DOA
    % Estimation// book)
    % A is the APM matrix with rows for each element (Mxd), d is signals? (pg 38)
    
    % Here A is A(theta), and P is a function of theta (thus loop over theta?)
    % do this in a subfunction to keep notation simple?
    % Note: ' is the Hermitian transpose, which is intended
    P = A*(inv(A'*A))*A'; % can I remove the inverse?
        
    tr(n) = trace(P*C);
        
    
end



% Ok ... what have I got?
[mxtr,ix] = max(tr); 

end

function [mxtr,ix] = find_mle_dual(APM,C,ii)
% FIND MLE - loop over APM bearings, finding the best fit to the signal

nt = length(APM.BEAR);

% TWO BEARING SEARCH
% create the APM excluding the ix bearing
B  = subsref_struct(APM,setdiff(1:nt,ii),nt,1); % < --- SLOW ... could get rid of this, use indexing instead

% get the shortened index length
nt = length(B.BEAR);

tr = NaN(nt,1);

    
for n = 1:nt
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    % make it a column vector (M element x d signals), M = 3
    % A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
    
    % CONSTRUCT THE APM MATRIX
    % Construct the two signal version
    A = [ APM.A13R(ii)+ 1i.*APM.A13I(ii)  APM.A23R(ii) + 1i.*APM.A23I(ii)  1+1i.*0; ...
            B.A13R(n) + 1i.*B.A13I(n)       B.A23R(n)  + 1i.*B.A23I(n)     1+1i.*0; ].'; % < --- SLOW
        
    % keyboard    
    
    % max(trace(P*C))
    % C is the covariance matrix (R in the //Introduction to DOA
    % Estimation// book)
    % A is the APM matrix with rows for each element (Mxd), d is signals? (pg 38)
    
    % Here A is A(theta), and P is a function of theta (thus loop over theta?)
    % do this in a subfunction to keep notation simple?
    % Note: ' is the Hermitian transpose, which is intended
    % This is eqn 11.b in Ziskind and Wax 1988
    % P = A*(inv(A'*A)*A'); % can I remove the inverse?   % < --- SLOW
    P = A*( (A'*A)\A' ); % inverse removed
    
    tr(n) = trace(P*C); % eqn 13.a  % < --- SLOW, maybe a faster way is tto just compute diag elements by hand?
        
    
end



% Ok ... what have I got?
[mxtr,iy] = max(tr); 


% % TEST PLOT
% figure
% plot(B.BEAR, real(tr)), hold on,
% plot(B.BEAR(iy),real(tr(iy)),'r*')


% ... so my two solutions so far are ix, and ii
% need to translate ix into the APM index

% TRANSLATE IX TO ABS INDEX
ix = find(APM.BEAR == B.BEAR(iy));


% TO DO
% look at A, P at the two bearings with really high tr and make sure the
% math is working right ...

end

function  [absix,TR] = mle(APM,C)
% MLE ALGORITHM - DOA via brute force search 
% absix = mle(APM,C)
%
% INPUTS
% APM  - Standard APM struct
% C    - data covariance matrix (from cross spectra struct)
% opt  - (OPTIONAL) plot the 2 bearing search surface
%
% OUTPUT
% absix - absolute index of the APM [single dual dual] at the estimated
%         direction of arrival
% TR    - matrix of trace values for the 2 bearing search map
%

% TO DO
% What should this output?
% Test for single/dual? 
% 
% 
% probably efficiencies to be had pre-computing all the P matricies first
% and storing them in the APM struct, maybe the 3d version for the two
% signal case?
%
% Other improvements based on refs that cite this method?
% More extensive study?: why would singular C matrix cause problems for
% Don? Look at 1983 ref?
%
% Parse this out into functions, move music.m into hf_doa_methods?
%
% MLE/MUSIC outline:
% load/get CS
% get first order region (whole CS or range cell?)
% run doa (loop over range cells?)
% run single dual tests
% make the radial file 

% tic

% check for test case
if strcmp('--t',APM), test_case; return, end

% SINGLE BEARING SEARCH
% Find the single bearing solution
% will probably need the value of the trace to compare with dual brg
% solution
[tr,ib] = find_mle(APM,C);



% TWO BEARING SEARCH
% now, combine this (one brg) soln with the other elements of the apm,
% compute and store the trace
n = length(APM.BEAR);

TR = NaN(n,n);

for ix = 1:n
   
    % just brute force compute the vaule of the trace for the whole grid of
    % possbile two bearing solutions
    %
    % This is slightly modified from the simlarly named subfunction in
    % mle_ap.m
    TR(:,ix) = find_mle_dual(APM,C,ix);
    
    
end

[~,r,c] = max_2d(TR);

absix = [ib r c];




end


function [mxtr,ix] = find_mle(APM,C)
% FIND MLE - single bearing case
%


% Get number of bearings to loop over
nt = length(APM.BEAR);

tr = NaN(nt,1);



% As a column vector for each n
% not really faster to do this here
A = [ APM.A13R(:)+ 1i.*APM.A13I(:)  APM.A23R(:) + 1i.*APM.A23I(:)  1+1i.*zeros(size(APM.A13R(:))) ].';

for n = 1:nt
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    % make it a column vector (M element x d signals), M = 3
    % A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
    
%     % PRE CONSTRUCT THE APM MATRIX
%     % If nargin == 2 we're looking for a 2 brg solution and need to modify the
%     % APM search    
%     A = [ APM.A13R(n)+ 1i.*APM.A13I(n)  APM.A23R(n) + 1i.*APM.A23I(n)  1+1i.*0 ].';

    
    
    % max(trace(P*C))
    % C is the covariance matrix (R in the //Introduction to DOA
    % Estimation// book)
    % A is the APM matrix with rows for each element (Mxd), d is signals? (pg 38)
    
    % Here A is A(theta), and P is a function of theta (thus loop over theta?)
    % do this in a subfunction to keep notation simple?
    % Note: ' is the Hermitian transpose, which is intended
%    P = A*(inv(A'*A))*A'; % can I remove the inverse?
    P = A(:,n)*(inv(A(:,n)'*A(:,n)))*A(:,n)';

    tr(n) = trace(P*C);
        
    
end



% Ok ... what have I got?
[mxtr,ix] = max(tr); 

end

function tr = find_mle_dual_off_by_one(APM,C,ii)
% FIND MLE - loop over APM bearings, finding the best fit to the signal
%
% BRUTE FORCE VERSION

% BRUTE FORCE VERSION

% nt = length(APM.BEAR);
% 
% TWO BEARING SEARCH
% create the APM excluding the ix bearing
% B  = subsref_struct(APM,setdiff(1:nt,ii),nt,1);

% % get the shortened index length
% nt = length(B.BEAR);
nt = length(APM.BEAR);

tr = NaN(nt,1);

    
for n = 1:nt
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    % make it a column vector (M element x d signals), M = 3
    % A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
    
    % CONSTRUCT THE APM MATRIX
    % Construct the two signal version
    A = [ APM.A13R(ii)+ 1i.*APM.A13I(ii)  APM.A23R(ii) + 1i.*APM.A23I(ii)  1+1i.*0; ...
          APM.A13R(n) + 1i.*APM.A13I(n)   APM.A23R(n)  + 1i.*APM.A23I(n)   1+1i.*0; ].';
        
        
    
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

keyboard


return

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

function tr = find_mle_dual(APM,C,ii)
% FIND MLE - loop over APM bearings, finding the best fit to the signal

nt = length(APM.BEAR);

% TWO BEARING SEARCH
% create the APM excluding the ix bearing
B  = subsref_struct(APM,setdiff(1:nt,ii),nt,1);

% get the shortened index length
nt = length(B.BEAR);

tr = NaN(1,nt);

    
for n = 1:nt
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    % make it a column vector (M element x d signals), M = 3
    % A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
    
    % CONSTRUCT THE APM MATRIX
    % Construct the two signal version (#antennas x #sources)
    A = [ APM.A13R(ii)+ 1i.*APM.A13I(ii)  APM.A23R(ii) + 1i.*APM.A23I(ii)  1+1i.*0; ...
            B.A13R(n) + 1i.*B.A13I(n)       B.A23R(n)  + 1i.*B.A23I(n)     1+1i.*0; ].';
        
        
    
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



% OUTPUT trace vector with same indexing as APM
tr = [tr(1:ii-1) NaN tr(ii:end)]';


return


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



% DEV TEST
function test_case

[C,APM] = test_data_tonys;

tic
% MLE FULL SEARCH
[ix,TR] = mle(APM,C);
toc

tic
% MLE AP FOR COMPARISON
[absix,AP] = mle_ap(APM,C);
toc

if isequal(ix(:),absix(:))
    disp('mle vs mle_ap ... get same result!')
else
    disp('mle vs mle_ap ... INCONSISTENT!')  
end

keyboard

% PLOT THE 2 BRG SEARCH SURFACE

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




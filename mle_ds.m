function  [absix,TR,delta] = mle_ds(APM,C)
% MLE DS - MLE for Distributed source via alternating projection
% absix = mle_ds(APM,C)
%
% REFERENCE
% An amalgam of Ziskind and Wax, 1988, Van Trees 2002, Meng et al. 1996 and
% Read, 1999
%
%
% INPUTS
% APM  - Standard APM struct
% C    - data covariance matrix (from cross spectra struct)
%
% OUTPUT
% absix - absolute index of the APM [single dual dual] at the estimated
%         direction of arrival

% Copyright (C) 2016 Brian M Emery
%
% Version Sept 2016
%
% From mle_ap.m

% TO DO
% - Need to handle the case of bearings near the edges
% - Combine different model-data covariances for more signal soures
% - Store the model data with the APM (and delta array)
% - Will likely want to do alternating projection for single bearing case


% check for test case
if strcmp('--t',APM), test_case;  end


% number of bearings
n = length(APM.BEAR);

% Construct the matrix of antenna location vectors 
% length(APM.BEAR) columns of px1 vectors, p = number of elements
A = [ APM.A13R(:)+1i*APM.A13I(:)  APM.A23R(:)+1i*APM.A23I(:)  ones(n,1)+1i*zeros(n,1) ].';

% Form the distributed source model-data covariance matricies ...
% for every value of delta, but for a single signal source
[Cmd,delta] = compute_model_cov(A,APM.BEAR);


% Now Combine different model-data covariances for more signal soures
% ... see Read, 1999 equation 31 


% SINGLE BEARING SEARCH
% Find the single bearing solution
% will probably need the value of the trace to compare with dual brg
% solution
[tr,ix] = find_mle(A,Cmd,C);


keyboard


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


if length(absix > 3)
    
    % grab the three (1 single and 2 dual) bearing solutions
    absix = absix([1 end-1 end]);
else
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


function [mxtr,ix] = find_mle(A,S,C)
% FIND MLE - single bearing case
%
% 


% number of bearings and deltas
[nbrg,ndel] = size(S);

% container for liklihood function
L = NaN(nbrg,ndel);


% loop over bearings
for b = 1:nbrg
    for d = 1:ndel
        
        if ~isempty(S{b,d})
            % Compute the liklihood function
            % ... See Van Trees, 2002 p 986 (eqn 8.280) and Read, 1999 eqn 10.
            % They have similar elements, Read has trace(R*inv(C)) while Van Trees
            % has trace(inv(C)*R) ...
            L(b,d) = - ( log(det(S{b,d})) + trace(S{b,d}\C) );
            
            
            %     % max(trace(P*C))
            %     % C is the covariance matrix (R in the //Introduction to DOA
            %     % Estimation// book)
            %     % A is the APM matrix with rows for each element (Mxd), d is signals? (pg 38)
            %
            %     % Here A is A(theta), and P is a function of theta (thus loop over theta?)
            %     % do this in a subfunction to keep notation simple?
            %     % Note: ' is the Hermitian transpose, which is intended
            %     P = A*(inv(A'*A))*A'; % can I remove the inverse?
            %
            %     tr(n) = trace(P*C);
            %
        end
    end
end

keyboard

% Ok ... what have I got?
[mxtr,ix] = max(tr); 

end

function [mxtr,ix] = find_mle_dual(APM,C,ii)
% FIND MLE - loop over APM bearings, finding the best fit to the signal

nt = length(APM.BEAR);

% TWO BEARING SEARCH
% create the APM excluding the ix bearing
B  = subsref_struct(APM,setdiff(1:nt,ii),nt,1);

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
            B.A13R(n) + 1i.*B.A13I(n)       B.A23R(n)  + 1i.*B.A23I(n)     1+1i.*0; ].';
        
        
    
    % max(trace(P*C))
    % C is the covariance matrix (R in the //Introduction to DOA
    % Estimation// book)
    % A is the APM matrix with rows for each element (Mxd), d is signals? (pg 38)
    
    % Here A is A(theta), and P is a function of theta (thus loop over theta?)
    % do this in a subfunction to keep notation simple?
    % Note: ' is the Hermitian transpose, which is intended
    % This is eqn 11.b in Ziskind and Wax 1988
    P = A*(inv(A'*A))*A'; % can I remove the inverse?
        
    tr(n) = trace(P*C); % eqn 13.a
        
    
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

% construction bits

function [H,delta] = compute_H(A,th)
% COMPUTE H - function of phi (theta and delta) 
% Pre compute H (equation 53)
% 
% INPUT
% A(theta) - the matrix of px1 location vector of antenna array (p = #elements) 
%            ... for each theta
% th       -  (aka theta) range of theta given by APM

disp('computng H ...')

% DEFINITIONS

% the range of distributed source widths
delta = 0:40;

% number of bearings
n = length(th);

% number of elements. I use p for this, and g for the ds shape function
p = size(A,1);


% Preallocate H 
% Starting with the single bearing case, it appears that H is going to be
% size n by (length delta) ie the number of widths 
% delta = 0:40, theta given by the APM
% 
% Note that each column of this is like the MUSIC DOA function, in fact the
% first column (delta = 0) should be the same as the MUSIC DOA
% H = zeros(n,length(delta));
H = cell(n,length(delta));

% define container for insides of integral
I = zeros(p,p,n);

% Now populate H, computing the 'integral' of eqn 53
%
% Two outer loops define the value of psi for which we comptue H

for b = 5:n-5                   % loop over bearings 
    for d = 1:length(delta)     % loop over deltas
        
        
        % COMPUTE H FOR EACH PSI
        
        % Get source distribution shape. This gets applied to each bearing
        % ... and its a function of theta
        % not change in notation from the paper, which has 2 definitions of
        % p
        g = src_dist(th,th(b),delta(d)); 
        
        if d == 7, keyboard, end
        
        % Individually compute the products inside the integral, then add
        % them all up
        for t = 1:n
            
            % elements in side integral of Equation 53 
            I(:,:,t) = A(:,t)*g(t)*A(:,t)';
            
        end
        
        
        % this needs indexing of some kind on H *** ONLY COMPUTING 1 H!!
        % **
        H{b,d} = sum(I,3);
        
        % zero out container
        I = zeros(p,p,n);
                
    end
end

disp('... done')

end



function p = src_dist(theta,theta_i,delta)
% SOURCE DISTRIBUTION - equation 68 in Valaee et al, 1995; butterworth form
%
% INPUTS
% theta     - bearing (free paramter) 
% theta_i   - source centroid location
% delta     - source 3db width
%
% Output is unitless, as long as input units (e.g. degrees) are consistent
%
% % equn 68 for example source distribution
% p = @(theta)(1./(1+((theta-thetai)./5).^2));
% x = -20:20; plot(x,p(x));

% see also gaussian_source.m


K = 1;

p = (K.^2) ./ ( 1+ ( (theta-theta_i)./delta ).^2 );


% if delta == 0, p is a delta function .. I think
if delta == 0
    p(isnan(p)) = (K.^2); 
end

end


function test_case
% TEST CASE
% 

% From Valaee et al. 1995, section V
% 
% CD? sources
% 20 element PA, theta = [10 13], delta = [1 2]


% ... simpler test: just get one of the model data and use that as data?
% SIMLUATION SOME DISTRIB SOURCE DATA

% Create the idealized pattern 
APM = make_ideal_pattern(225, 1:360);

%S{180,20}
C =[ 0.5000    0.3919    0.6653
    0.3919    0.5000    0.6653
    0.6653    0.6653    1.0000];

[ix,delta,DOA] = mle_ds(APM,C);


keyboard


% Specify a point source at bearing 205 
R.RadComp =  300 ;
R.BEAR =  205 ;


% Get defaults
CFG = radar_sim_config;

% override some defaults
CFG.wind_noise = 0;    
CFG.sigl_noise = 1;    
CFG.curr_noise = 0;    
CFG.SNR = 20;


[CSA, ~] = radar_simulation_ds(APM,CFG,R.BEAR,R.RadComp);

% Get freq bin by hand 
fx_ds = 422;
Cds = make_cov(CSA,fx_ds,1); 

% NOW TEST MLE
[ix,delta,DOA] = mle_ds(APM,Cds);


keyboard



h = surf(delta,APM.BEAR,real(10*log10(DOA))); set(h,'edgecolor','interp');

keyboard


% TEST CASE - TONYS DATA

% Dual bearing test case, ideal APM
% [C,APM] = tonys_test_data;
%
% Create the covariance matrix from de Paolo's example:
C=[ 0.2162          0.0303-0.0090i  0.3170-0.0063i; ...
    0.0303+0.0090i  0.0436         -0.0091+0.0213i; ...
    0.3170+0.0063i -0.0091-0.0213i  0.5416];

% Create the idealized pattern ... 1 deg here 
APM = make_ideal_pattern(225, 1:360);

[ix,delta,DOA] = dspe(APM,C);


% this appears to compare well with MUSIC! (205 and 330) 

h = surf(real(10*log10(DOA))); set(h,'edgecolor','interp');

plot3([0 0; 0 0 ], [205 330; 205 330], [0 0; 40 40],'-r')





end


function [H,delta] = compute_model_cov(A,th)
% COMPUTE MODEL COVARIANCE - function theta and source width, delta
%
% Compute the covariance of the 'model data' for a distributed source
% 
% Following Read, 1999 equation 27
% 
% INPUT
% A(theta) - the matrix of px1 location vector of antenna array (p = #elements) 
%            ... for each theta
% th       -  (aka theta) range of theta given by APM
%
% OUTPUT
% H        - In Read, 1999, this is Cmm which is the 

disp('Computing model-data ...')

% DEFINITIONS

% the range of distributed source widths
delta = 1:40;

% number of bearings
n = length(th);

% number of elements. I use p for this, and g for the ds shape function
p = size(A,1);


% Preallocate H 
% Starting with the single bearing case, it appears that H is going to be
% size n by (length delta) ie the number of widths 
% delta = 0:40, theta given by the APM
% 
% Note that each column of this is like the MUSIC DOA function, in fact the
% first column (delta = 0) should be the same as the MUSIC DOA
% H = zeros(n,length(delta));
H = cell(n,length(delta));

% define container for insides of integral
C = zeros(p,p,n);


% Now populate H, computing the 'sum' of eqn 27
%
% Two outer loops define the value of psi for which we comptue H

for b = 5:n-5                   % loop over bearings (source locations)
    for d = 1:length(delta)     % loop over deltas
        
        
        % COMPUTE Cmm for each valute of delta and theta
        
        % Get source distribution shape. This gets applied to each bearing
        % ... and its a function of theta
        % not change in notation from the paper, which has 2 definitions of
        % p
        g = gaussian_source(th,th(b),delta(d));
                
        % to figure out: sum this over all bearings, but need loop index
        % for the source bearing location ... 
        
        % sum performed over all the bearings (aka the grid), "sum of the
        % individual covariances of the component point-source signals"
        % most of these will be zero ... so a possible speed up is to just
        % find the relevant ones to compute and sum over ...
        for ii = 1:n
            C(:,:,ii) =  g(ii).'.* A(:,ii)*A(:,ii)';
        end
        
        
        % When you add these all up, you get the cov for one distributed
        % source, with a given delta, at one center bearing ...
        H{b,d} = sum(C,3);
        
        
        % zero out container
        C = zeros(p,p,n);
                
    end
end

disp('... done')

end


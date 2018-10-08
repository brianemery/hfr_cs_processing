function err = music_error(A,C,th,K,n,ix) %,SNR)
% MUSIC ERROR - standard deviation of DOA error from Stoica and Nehorai 1989
% err = music_error(A,C);
% 
% Gives the error variance of a MUSIC solution based on equation 3.12 of 
% MUSIC, Maximum Likelihood and the Cramer-Rao Bound, by Stoica and 
% Nehorai, 1989
% 
% INPUT
% A   - the complex array matrix (e.g. from get_array_matrix.m)
% C   - the data covariance matrix
% th  - bearing (in degrees) associated with A
% K   - number of data snapshots
% n   - number of signals assumed
% ix  - the index of the DOA estimate
% % SNR - (dB) external SNR estimate to use instead of eigenvalues for SNR ...
%
% OUPUT
% err - Error standard deviation in degrees
% 
% EXAMPLE
%
% NOTES
% IF:
% mse    - Mean square error of MUSIC direction estimate (see below) 
% evar   - Variance of MUSIC direction estimate (from Stoica and Nehorai,
%          1989)
% ebias  - Bias of MUSIC direction estimate (from Xu and Buckley 1992)
%
% Then  MSE = evar + (bias)^2


% Copyright (C) 2016 Brian Emery
%
% Version 29-Nov-2016 15:48:16

% TO DO
% - Fix for ULA test:
%   -try using average of noise eigenvectors for sigma
%   -build test using ULA
%
% - This holds for sufficiently large K (N in the paper). How few is too
%   few? Need to test (K > m antennas in paper)
%
% DONE
% - Investigate the use of the smallest eigen value as an estimate of the
%   noise level (See NOTES in test_case) 10*log10(D(3)/D(1))
% - implement 7.5a instead and compare (not really possible - theory only)
% - Try second part of 3.12 instead?
% - better estimate of sigma based on SNR doesnt seem to work, c.f.
%   plot_two_brg_music_error.m
% - Ouput an error for each of the ix bearings?
% - paper says 'unit norm' eigenvectors - eig outputs these by default


% check for test case
if strcmp('--t',A), test_case, return, end


% GET EIGEN VECTORS AND EIGEN VALUES
% note: "[V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
% full matrix V whose columns are the corresponding eigenvectors so
% that X*V = V*D.", ie C*eigVectors=eigVectors*eigValues
%
%  [...] = eig(...,'vector') returns eigenvalues in a column vector 
%  instead of a diagonal matrix.
[V,D] = eig(C); %,'vector'); 

% sort largest to smallest to be like the paper
[D,edx] = sort(diag(D),'descend');
V = V(:,edx); clear edx

% Get number of columns in covariance matrix
m = size(C,2);


% COMPUTE da/dw
% Array derivative at each value of theta
da = diff_centered(A,th*pi/180);


% FORMULATE EQUATION 3.12  

% compute G (matrix of noise eigen vectors)
G = V(:,(n+1):m);


% Compute noise power 
% Note that equation 3.2 has the noise level
% (sigma) given by the noise eigen values. So I will use that here
%
% Mean of noise eigenvalues suggested by Krim and Viber 1996
sigma =  mean(D(n+1:end)); % min(D); % or? (non-dB)   

% compute U (see equation 3.9). 
U = zeros(m);

% loop over signals
for i = 1:n
    U = U  + ( D(i)/((sigma-D(i)).^2) ) * ( V(:,i)*V(:,i)' );
end

U = sigma.*U;


% Run this part of the calc at each DOA
err = NaN(size(ix));

for i = 1:numel(ix)
    
    % compute H ( h(w) in the paper ) 
    % Big err results when this is small, derivative near zero? 
    H = da(:,ix(i))'* (G*G') * da(:,ix(i)); % <----- CHECK FOR MISTAKE, should this be G_k? not whole matrix?
    
    % Equation 3.12
    err(i) = (1/(2*K)).* (( A(:,ix(i))' * U * A(:,ix(i))) ./ H); 
    
    % keyboard
    % Equation 3.12, second part (would go here?)
        
end

% I suspect that the output is a variance in radian squared, so convert
% this to a standard deviation in degrees:
err = rad2_to_deg(real(err));

% if err > 100,
%     
%     % look at some things (is it APM or something else?)
%     for j = 1:numel(th), DA(j) = da(:,j)'* da(:,j); end
%     plot(th,DA,'-b.'), hold on
%     plot(th(ix),DA(ix),'r*')
%     
% end

end


function y = db2pow(ydB)
%DB2POW   dB to Power conversion
%   Y = DB2POW(YDB) converts dB to its corresponding power value such that
%   10*log10(Y)=YDB
%
%   % Example:
%   %   Convert 12dB to Power.
%
%   y = db2pow(12)      
%

%   Copyright 2006 The MathWorks, Inc.


y = 10.^(ydB/10);

end

function test_case
% TEST CASE
% 
% TO DO 
% run some simulations repeating the figures in MUSIC, Maximum Likelihood, 
% and Cram&-Rao Bound: Further Results and Comparisons 
% see Section VII for details
%
% 
% Makes a figure showing the DOA function vs bearing just like fig 9 in the
% De Paolo and Terrill Scripps report
%
% Lots of code from music.m
%
% NOTES
% in check_two_brg_understanding.m I looked at the modelled SNR vs SNR
% determined by the ratio of eigen values, using this:
% [V,D] = eig(C); D = diag(D);
% SNR(j)
% 10*log10(D(3)/D(1))
% w
% snr_123 = get_SNR(CSA,fbin)
% ... and I get an answer that is quite close to the modelled, particularly
% given the different ways of computing SNR

% (APM,DOA,singleIdx,dualIdx)

% % MAKE SOME DATA 
% % This is the 'bad data' case - see check_two_brg_understanding
% 
% SNR = 10;
% brgs = [205 235]; % 240];
% 
% A = make_ideal_pattern(249,90:270);
% 
% % FROM RADAR SIMULATION
% 
% % standard stuff ...
% CFG = configure_radar;
% CFG = configure_sim(CFG);
% 
% 





% USE TONYS DATA

% Create the covariance matrix from de Paolo's example:
% NOTE that he made this from a simulation with currents input at 205 and
% 330 degrees which MUSIC gets kind of wrong
C=[ 0.2162          0.0303-0.0090i  0.3170-0.0063i; ...
    0.0303+0.0090i  0.0436         -0.0091+0.0213i; ...
    0.3170+0.0063i -0.0091-0.0213i  0.5416];

th =  1:360;
K = 3;
n = 2;

% Create the idealized pattern 
APM = make_ideal_pattern(225,th);

% Get array matrix
A = get_array_matrix(APM);

% -- put music.m here ---

[DOA,ix] = music(A,C,th,2);

% % GET EIGEN VECTORS AND EIGEN VALUES
% % note: "[V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
% % full matrix V whose columns are the corresponding eigenvectors so
% % that X*V = V*D.", ie C*eigVectors=eigVectors*eigValues
% [eigVectors,eigValues]=eig(C);
% 
% % Compute DOA's
% [DOA,singleIdx,dualIdx] = getDOAs(APM,eigVectors,eigValues);
% 
% -- -- -- -- -- -- -- -- 

% Run the calculation
% err = music_error(A,C,th,K,n,ix) A,C,th,K,n,ix,SNR
err1 = music_error(A,C,th,K,1,ix{1});

% Now make sure two bearing case works
err2 = music_error(A,C,th,K,n,ix{2});

% cross check
[err_var,~] = music_error2(A,C,th,K,n,ix{2});

rad2_to_deg(err_var)
 keyboard

%  % interesting ... but meaningless?
%  for ii = 1:numel(th), err1(ii) = music_error(A,C,th,K,1,ii); end
%  figure, plot(th,err1)
 


% OLDER BELOW - USES (probabaly incorrect) EXTERNAL SNR EST METHOD

% Compare with CRB
A = get_array_matrix(APM);

SNR = 2:35;

for i = 1:length(SNR)
    % guess what the SNR might be
    err1(i) = music_error(A,C,APM.BEAR,3,1,singleIdx,SNR(i));
    
    err2(i) = music_error(A,C,APM.BEAR,3,1,singleIdx,SNR(i));
    
end

CRB2 = compute_crb(A,APM.BEAR,SNR,dualIdx,1);

CRB1 = compute_crb(A,APM.BEAR,SNR,singleIdx,1);

keyboard
figure, hold on
h1 = plot(SNR,CRB1,'r');
h2 = plot(SNR,CRB2,'b');
h3 = plot(SNR,err1,'-go');
h4 = plot(SNR,err2,'-m*');

legend([h1 h2 h3 h4],'CRB n=1','CRB n=2','\sigma_M_U_S_I_C n=1','\sigma_M_U_S_I_C n=2')



end


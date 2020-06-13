function [err_var,bias] = music_error2(A,C,th,K,n,ix)
% MUSIC ERROR - RMS DOA error (bias and variance) for MUSIC
% [err_var,bias] = music_error2(A,C,th,K,n,ix)
% 
% Computes the error variance of a MUSIC solution based on equation 3.12 of 
% 'MUSIC, Maximum Likelihood and the Cramer-Rao Bound', by Stoica and 
% Nehorai, 1989, and the error bias following 'Bias Analysis of the MUSIC
% Location Estimator' by Xu and Buckley, 1992
% 
% INPUT
% A   - the complex array manifold (e.g. whole thing from get_array_matrix.m)
%       rows = numer of elements, columns length(th)
% C   - the data covariance matrix
% th  - bearing (in degrees) associated with A
% K   - number of data snapshots
% n   - number of signals assumed
% ix  - the index of the DOA estimate (index of th)
%
% OUPUT
% err_var - Error variance in radians^2
% bias    - MUSIC bias error in radians(?)
%
% Then ...
% MSE = err_var + bias^2
%
% STATUS
% The test case to this function produces comparable figures to XB92
% figures, providing confirmation that both the bias calculation and
% variance calculation are correct (or as correct as XB92). 
% 
% ... but I dont reproduce the XB92 figure 4. Look carefully at how I'm
% making the plot, also now taking mean of noise eigenvalues but this
% doesnt seem to help


% Copyright (C) 2017 Brian Emery
%
% Version 09-Feb-2017 11:37:47
% from music_error.m

% TO DO
% - Still problems reproducing figures. Found a mistake with the summation
%   (SM below) but this seems to make little difference. Need to do more
%   vetting ...
%
% - This holds for sufficiently large K (N in SN89). How few is too
%   few? Need to test (K > m antennas in paper)
% - Check that in the bias calculation the a(theta)/norm(a(theta),2) is
%   equivalent to what I have for A
%
% DONE
% - Investigate the use of the smallest eigen value as an estimate of the
%   noise level (See NOTES in test_case) 10*log10(D(3)/D(1))

% DEFINITIONS
% A   - the complex array matrix (e.g. from get_array_matrix.m)
% C   - the data covariance matrix
% th  - bearing (in degrees) associated with A
% K   - number of data snapshots
% n   - number of signals 
% ix  - the index (of th) of the DOA estimate
% V   - matrix of eigenvectors
% D   - matrix of eigenvalues


% check for test case
if strcmp('--t',A), test_case, return, end


% GET EIGEN VECTORS AND EIGEN VALUES

% Get eigen-decomposition
[V,D] = eig(C); 

% sort largest to smallest to be like the paper
[D,edx] = sort(diag(D),'descend');  
V = V(:,edx); 


% % note: "[V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
% % full matrix V whose columns are the corresponding eigenvectors so
% % that X*V = V*D.", ie C*eigVectors=eigVectors*eigValues
% %
% %  [...] = eig(...,'vector') returns eigenvalues in a column vector 
% %  instead of a diagonal matrix.
% [V,D] = eig(C,'vector'); 
% 
% % sort largest to smallest to be like the paper
% [D,edx] = sort(D,'descend');
% V = V(:,edx); clear edx

% Get number of columns in covariance matrix (Array elements)
m = size(C,2);


% NORMALIZE A
% see just below eqn 4 - this is my interpretation, that is, I think a(th)
% is a vector in their definition of v(th)
for i = 1:length(th) 
    A(:,i) = A(:,i)./norm(A(:,i),2); 
end



% COMPUTE ARRAY DERIVATIVES

% Array derivative at each value of theta
da = diff_centered(A,th*pi/180);

% get array second derivative (v dot dot)
d2a = diff_centered2(A,th*pi/180);




% CALCULATIONS

% Eqn 3.12 SN89
% this is the error variance, same as original code in music_error.m
err_var = music_err_var(V,D,n,m,A,da,ix,K);


% Eqn 22 XB92
bias = music_bias(V,D,n,m,A,da,d2a,ix,K,err_var);

% 
% % I suspect that the output is a variance in radian squared, so convert
% % this to a standard deviation in degrees:
% err = rad2_to_deg(real(err));
% 

% estimate array snr

end

function bias = music_bias(V,D,n,m,A,da,d2a,ix,K,err_var)
% MUSIC BIAS - From Xu and Buckley, 1992
% 
% Paper notation: K sensors, D signals, v is normalized a(theta)
% N snapshots


bias = NaN(size(ix));

% Loop over each signal, compute bias for each (theta_d in paper)
for i = 1:numel(ix)
    
    bias(i) = eqn_22(V,D,n,m,A,da,d2a,ix(i),K,err_var(i));
    
        
end

% take real parts
bias = real(bias);

end

function bias = eqn_22(V,D,n,m,A,da,d2a,ix,K,err_var)
% EQN 22 - equation 22 of XB92, built in parts
% 
% Note change in inputs from XB92 notation:
% K sensors -> m
% D signals -> n
% v is normalized a(theta)
% vdot -> da
% N snapshots -> K
% e eigenvectors -> V
%
% Here, D is eigenvalues vector
% 

% PRELIMINARIES

% get matrix of noise eigenvectors En, and signal eigenvectors Es
%En = V(:,(n+1):m);
Es = V(:,1:n);

% get noise variance (sigma^2 in XB92, sigma in SN89)
% Mean of noise eigenvalues suggested by Krim and Viber 1996
% -- this needs to be min(D) to reproduce fig 3 in XB92
sigma2 =  min(D); %mean(D(n+1:end)); %  or? (non-dB)   

% compute 2nd derive of D_mu in paper -  eqivalent to I-Es*Es':
% (En*En') - (eye(m)-Es*Es')
% 
% ans =
% 
%    1.0e-15 *
ddD_mu = 2*da(:,ix)'* (eye(m)-Es*Es') * da(:,ix);


% ASSEMBLE THE PARTS


% build up the summation
% loop over assumed number of signals to build it up
SM = 0;
for k = 1:n   
        
    % hellow, K-D-1 = 0, unless 3-1-1 = 1, ... if dual then contrib from
    % this term
    % Check D is vector
    SM = SM + (( (m-n-1).*D(k).*sigma2 ) ./ (D(k) - sigma2).^2 ) * ...
                                 real( da(:,ix)'*V(:,k)*V(:,k)'*A(:,ix) );
                             
end



% Compute the first term
T1 = -(2/K)*SM/ddD_mu;


% compute the second term (see eqn 23)
% 6's cancel (top part in real( ) is D-triple-dot)
T2 = (real( d2a(:,ix)' * (eye(m)-Es*Es') * da(:,ix) )/ddD_mu ) * err_var;

 

% put it all together
% T1 has the minus sign already
bias = T1 - T2;

    
    
end

function err = music_err_var(V,D,n,m,A,da,ix,K)
% FORMULATE EQUATION 3.12 SN89
% this is the error variance
%
% compute it separately form the bias, and output what's needed

% compute G (matrix of noise eigen vectors)
G = V(:,(n+1):m);


% compute noise power 
% Note that equation 3.2 has the noise level
% (sigma) given by the noise eigen values. So I will use that here
% if nargin < 7
sigma = mean(D(n+1:end)); % min(D); % or ? (non-dB)
% else
%     % get noise level relative to max signal, using SNR input (in dB)
%     % 10*log10(D(1)/D(3))
%     % 10*log10(D(1)/sigma)
%     %
%     sigma = max(D) ./db2pow(SNR);
% end

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
    H = da(:,ix(i))'* (G*G') * da(:,ix(i)); 
    
    % Equation 3.12
    err(i) = (1/(2*K)).* (( A(:,ix(i))' * U * A(:,ix(i))) ./ H);
    
    % keyboard
    % Equation 3.12, second part (would go here?)
        
end

% take real part
err = real(err); % << maybe not yet??


end

function test_case
% TEST CASE
%
% See also music_error_test.m, which replicates Van Trees 2002 or something

% 
% % % SIMPLER TEST
% % % XB92 Section V A. ... single source at theta = 0:80, K = 10, N = 20, 
% % % SNR = 0 dB (that is, norm(x)^2/var(x) = 1)
% % 
% th = -90:90;
% M = 10;  % # elements in array
% K = 20;  % # snapshots
% n = 1;   % # signals
% SNR = 0; % dB
% 
% % Make the array use their definition (just past eqn 24) 
% [A,~] = make_ula_pattern(th,M,'xu');
% 
% bias = NaN(1,81); %size(brg));
% 
% % loop over source bearings
% for brg = 0:80
% 
%     [~,ix,~] = intersect(th,brg);
%     
%     [C,~] = radar_simulation_basic(A(:,ix),K,SNR);
%     
%     % ... weird ... do I need music?
%     [evar(brg+1),bias(brg+1)] = music_error2(A,C,th,K,n,ix); 
% 
%     
% end
% 
% % Apr 2018: put means of noise eigenvalues in here which breaks the figure
% % priviously I had:
% %
% % Fig 3. PRETTY CLOSE !!!
% plot(0:80,bias*180/pi), axis([0 80 0 0.2])
% title('Xu and Buckley 1992 Figure 3')
% 
% keyboard




% COMPLEX VALIDATION TEST -  XB92 Figure 4

% BASIC NARROWBAND SIMULATION
% 14 and 16 deg equipower uncorrelated sources, 500 trials per, vary SNR, 10 element
% ULA, N = 20, 
th = -90:0.01:90; %th = 5:0.01:25; %wtf did I do this?
M = 10;  % # elements in array
K = 20;  % # snapshots
n = 2;   % # signals

[~,ix,~] = intersect(th,[14 16]);

[A,~] = make_ula_pattern(th,M,'xu');

SNR = 15:50;
Nsim = 100; %500; % " 500 trials per SNR" % they used 800 trial


% NOTE: here is some boilerplate for monte-carlo radar simulation
% with minor mods. C.f. two_brg_simulation_basic.m

% Initialize storage
% keep only one of the bearings as they do in the paper
[R.Bear, R.bias, R.evar] = deal(NaN(Nsim*length(SNR),1)); % was ,n));
[R.asnr, R.SNR         ] = deal(NaN(Nsim*length(SNR),1));


% place holder
t = 1;

for mc = 1:Nsim
    for i = 1:numel(SNR)
        
        % SNR goes 20:50, plot shows bias and STD for 14 deg source
        [C,asnr] = radar_simulation_basic(A(:,ix),K,SNR(i));
        
        % get music doa solution
        [~,idx] = music(A,C,th,2);
        
        % get the errors
        [err_var,bias] = music_error2(A,C,th,K,n,idx{2});
 
        
        % saving only data associated with source at 14 deg
        % need to be careful with indexing here. tx is the index assoc with
        % th, j is the index associated with idx{2}
        % [tx,j] = min(idx{2});
        
        % allow for possibilty of less than 2 solutions
        s = length(idx{2}); % if s<2, keyboard, end
        
        % save data
        R.Bear(t,1:s) = th(idx{2}); %(j));
        R.bias(t,1:s) = bias;
        R.evar(t,1:s) = err_var;
        R.asnr(t,1)  = asnr;
        
        R.SNR(t) = SNR(i);
                
        t = t+1;
    end
end


% DATA BINNING
% ** plot only data associated with source at 14 deg ***

% get array of source true bearing
Bear = repmat(th(ix(1)),Nsim*length(SNR),1);

% % For the results that have one NaN (single bearing case), put the doa
% % result closer to the source (that is, so that the 19.5 doesn't get
% % counted as a doa for the 0 deg source)
% rx = find(isnan(sum(R.Bear,2)) & (R.Bear(:,1) > 10) );
% R.Bear(rx,2) = R.Bear(rx,1);
% R.Bear(rx,1) = NaN;


% these are straight averages at each SNR bin
[bias,~,~,~,~] = binData(R.bias(:,1), R.asnr(:), 1, SNR);
[evar,~,~,~,~] = binData(R.evar(:,1), R.asnr(:), 1, SNR);

% Here, use the binning to compute the mean part of RMS
% (RMS is difference, squared, mean-ed, sqrt) ...
[rmse,~,~,~,~] = binData( (R.Bear(:,1)-Bear(:)).^2 , R.asnr(:), 1, SNR);

% ... now take sqrt
rmse = sqrt(rmse);

% convert out of radians
bias = bias.*180./pi;
estd = rad2_to_deg(evar);

% bias in radians?, var in radians squared
figure
h1 = plot(SNR,bias,'b-'); hold on
h2 = plot(SNR,estd,'r--'); 
h3 = plot(SNR,rmse,'-g*');

title('Xu and Buckley 1992 Figure 4')

axis([20 50 0 0.5])

legend([h1 h2 h3],'Analytical Bias','Analytical \sigma','Sim \sigma')


keyboard


% Make XB92 Figure 1
% this figure suggests music is more accurate than they show. Perhaps I
% have something wrong with my simulation code or with music.m?
%
% This reproduces their bias estimate though

figure, 
plot([20 50],[14 14],'--'), hold on
plot([20 50],[16 16],'--')
axis([15 50 13.5 16.5])

[brg1,~,~,~,~]=binData(R.Bear(:,1),R.SNR, 1,SNR);
plot(SNR,brg1,'r.'),

[brg2,~,~,~,~]=binData(R.Bear(:,2),R.SNR, 1,SNR);
plot(SNR,brg1,'r.'),


[bias1,~,~,~,~]=binData(R.bias(:,1),R.SNR, 1,SNR);
bias1 = bias1.*180./pi;

plot(SNR,14+bias1,'b*'),

[bias2,~,~,~,~]=binData(R.bias(:,2),R.SNR, 1,SNR);
bias2 = bias2.*180./pi;

plot(SNR,16+bias2,'bo'),

title('Xu and Buckley 1992 Figure 1(b)')


keyboard

end














function CRB = compute_crb(A,th,SNR,n,K) % dm,
% COMPUTE CRB - Cramer-Rao Bound on DOA estimation 
% CRB = compute_crb(A,th,SNR,n,K)
%
% Computes the Cramer-Rao Lower Bound on the DOA error. That is:
% 
%           MSE(theta) = VAR(theta) >= CRB (TF09, eqn 1.23)
%
% See Tuncer and Friedlander, 2009, A.1.3 pg 38 for theoretical development
% and see Ziskind and Wax, 1989 for comparison figure
%
% NOTES
% The CRB is a variance (TF09 and Brandwood, 2009) given in radians^2, 
% so to get degrees, this converts to degrees^2 and take the square root.
%
% Ziskind and Wax define SNR as 10 * log(s^2/sigma^2) (signal and noise
% powers). Thus use db2pow to convert out of dB.
%
% INPUT
% A   - the FULL complex array matrix (needed to compute derivatives)
% th  - bearing (in degrees) associated with A
% SNR - the signal to noise ratio in dB (can be a vector - loops over this
%       by default
% n   - the INDEX of the assumed emitter (signal) bearings
% K   - the number of data snapshots
%
% OUTPUT
% CRB - the Cramer-Rao Lower Bound as a standard deviation (in degrees)
%
% EXAMPLE USAGE
% 
% % Get array matrix
% AA = get_array_matrix(A);
% 
% % Get index of APM bearings of emitter locations
% [~,ia,~] = intersect(A.BEAR,brgs);
% 
% % Snapshots
% K = 5;
% 
% CRB = compute_crb(AA,A.BEAR,2:35,ia,K);
% 
% h3 = plot(hx(1),2:35,CRB,'--g');
%  
% SEE ALSO
% this function derived from compute_crb_n_signals.m and
% compute_crb_phased_array.m which were both development cases. An
% alternative method is used in compute_crb_sn89.m, which has more test
% examples also.

% Copyright (C) 2016 Brian Emery
%
% Version 14-Oct-2016 16:00:06

% Check for test case
if strcmp('--t',A), test_case, return, end


% COMPUTE ADOT
% Array derivative at each value of theta
adot = diff_centered(A,th*pi/180);



% Convert SNR out of db (see Stoica and Nehorai, 1990 for example)
% for the CRB calculation ( y = 10.^(ydB/10) )
SNR = db2pow(SNR);

CRB = zeros(size(SNR));

for i = 1:length(SNR)
    
    % Just pass the parts of A and adot that are needed
    CRB(i) = crb_for_each_snr(A(:,n),SNR(i),K,adot(:,n));

end

% Output the CRB in degrees (as a standard deviation)
CRB = rad2_to_deg(CRB);

end

function CRB = crb_for_each_snr(A,SNR,K,adot)
% CRB FOR EACH SNR - CRB calc with scalar SNR input
%
% Differs from the subfunction in compute_crb_n_signals in the use of
% diff_centered.m to compute the derivatives

% A is number of elements x number of emitters
[M,n] = size(A);


% TF09 EQN 1.108
% Construct the array covariance matrix (needed for the Fisher Information 
% Matrix).  Assume equal SNR for each signal (for now). 
R = zeros(size(A,1));

for i = 1:n
    R = R + ( SNR .* A(:,i)*A(:,i)' );
end

% I think this is outside the sum
R = R + eye(size(A,1));


    
% TF09 EQN 1.110
% need this at each theta (need to preallocate ... each of these are a MxM
% matrix
dRdth = complex(zeros(M,M,n));

for i = 1:n
    dRdth(:,:,i) = SNR .* ( (adot(:,i) * A(:,i)') + (A(:,i) * adot(:,i)') );
end




% TF09 EQN 1.107
F = zeros(n);

for i = 1:n
    for j = 1:n
   
        % Compute elements of F
        % F(i,j) = trace( inv(R)*dRdth(:,:,i) * inv(R)*dRdth(:,:,j) );
        % .. the matlab way
        F(i,j) = trace( ( R\dRdth(:,:,i) ) * ( R\dRdth(:,:,j) ) );
    end
end

% TF09 EQN 1.109
% X^(-1) is equivalent to inv(X). Also, units are radians^2 
% This produces a matrix (I believe a covariance matrix). Thus, we want the
% cross terms? Also, Stoica and Nehorai take the real part?
CRB = (1/K).*inv(F);


% I found one example that indicated that the CRB is found on the diagonal.
% Perhaps changing the SNR of the different emitters would result in
% different values along the diagonal. Also, the cross terms are very
% small, FYI.
CRB = real(CRB(1,1));



end

function test_case
% TEST CASE FROM case_tf09_fig_1pt13 subfunction to compute_crb_n_signals
% CASE STUDY: TF09 FIGURE 1.13
%
% CRB vs SNR for two emitters
%
% This gets the right answer for Separation less than about 3 degrees
% Not sure why ... but will have to figure it out at some point
%
%
% NOTE 
% it's not clear in the text how many snapshots they used, but this gets
% the right answer
K = 2;  
SNR = 30; 

% elements in array
M = 8;

% Half the emitters separation in degrees, zero is normal to array
del = [0.4:0.2:2.6 3:30]/2;

% Define total range of bearings
th = -90:0.1:90;

% Define Array Matrix
% Construct the matrix of antenna location vectors
% length(APM.BEAR) columns of px1 vectors, p = number of elements
[A,~] = make_ula_pattern(th,M); %,'krim');


% Run calc for each DEL by Looping over DEL
CRB = zeros(size(del));

for i = 1:length(del)
    
    % Define emitter bearing indecies (round to tenth of a degree)
    [~,ia,~] = intersect(round(th*10),round([-del(i) del(i)].*10));
    
    if isempty(ia), keyboard, end
    
    % Compute CRB
    CRB(i) = compute_crb(A,th,SNR,ia,K);
    
    
end



% figure 1.13
h1 = plot(2*del,CRB,'-m*'); hold on;
axis([0 30 0 0.7])
xlabel('Separation (deg)')
ylabel('Standard Deviation (deg)')
grid
title('TF 09 Figure 1.13')





% PLOT VS SNR

SNR = 2:2:30;

% For 8 element array
CRB8 = compute_crb(A,th,SNR,ia,K);

% SIMILARLY FOR SEASONDE

% RUN IT FOR SEASONDE TOO

% th is 10th of a degree so get points about 40 deg apart
ia = [1000 1400];

APM = make_ideal_pattern(225,th);

A = get_array_matrix(APM);

% compare with 3 element monopole
[A3,~] = make_ula_pattern(th,3); %,'krim');


% CRB CALCS

% Compute CRB
CRBSS = compute_crb(A,th,SNR,ia,K);

% single brg
CRBSS_1 = compute_crb(A,th,SNR,ia(1),K);

% 3 element ULA
CRB3 = compute_crb(A3,th,SNR,ia,K);


% MAKE THE FIGURE
figure(2)
h1 = plot(SNR,CRB8,'-b*'); hold on
h2 = plot(SNR,CRBSS,'-rs'); 
h3 = plot(SNR,CRBSS_1,'-y*'); 
h4 = plot(SNR,CRB3,'-ko');



hg = legend([h1 h2 h3 h4],'CRB (8 elem ULA)','CRB SS', 'CRB SS (n=1)', ...
              'CRB (3 elem ULA)');
          
          set(hg,'box','off')

          axis([0 35 0 35])
          
          
xlabel('SNR (dB)')
ylabel('Standard Deviation (deg)')
grid
title('CRB from TF09')

keyboard


% EXTRA CHECK

% see compute_crb_phased_array.m - single bearing case?

% Compute dbar^2 (1.69)
dbar2 =  sum(dm.^2) ;

SNRv = db2pow(SNR);

CRB_25 = 1./( 8 .* pi.^2 .* K .* SNRv .* cosd(th(ix(1))).^2 .* dbar2  );

h2 = semilogy(SNR, rad2_to_deg(CRB_25) ,'ro');




end


function pw = db2pow(db)

pw = 10.^(db/10);

end



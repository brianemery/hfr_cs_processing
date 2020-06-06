function [LR,gm,Rm,LR2] = glrt(A,R,K) 
% GLRT - signal detection with the generalized liklihood ratio test  
% [LR,gm,Rm,LR2] = glrt(A,R,K)
%
% Obtains the likely number of sources. Assumes K > M (snapshots >
% antennas). See Abramovich and Spencet 2004 for details).   
% 
% INPUTS
% A    - Array matrix at theta(s) given by the DOA solution
% R    - Data Covariance matrix
% K    - Number of data snapshots 
%
% OUTPUTS
% LR   - the Normalized Likelihood Ratio (method of Abramovich et al. 2006)
%        or the sphericity LR
% GM   - 
% Rm   - Model cov created from the DOA solutions
% LR2  - LR from Ottersten et al. 1993
%
%  ... (ideally) m   - the likely number of sources
% 
%
% REFERENCE
% Abramovich et al. 2006, Detection-Estimation of Gaussian Sources by
%   Nonuniform Linear Antenna Arrays: An Overview of the Generalized
%   Likelihood-Ratio Test Approach. Ch 2 of Advances in Direction of Arrival
%   Estimation, Chandran Editor (a book, from the library)
%
% @inproceedings{abramovich2004performance,
%   title={Performance breakdown of subspace-based methods in arbitrary antenna arrays: GLRT-based prediction and cure},
%   author={Abramovich, YL and Spencer, NK},
%   booktitle={Acoustics, Speech, and Signal Processing, 2004. Proceedings.(ICASSP'04). IEEE International Conference on},
%   volume={2},
%   pages={ii--117},
%   year={2004},
%   organization={IEEE}
% }
%
% Ottersten et al. 1993, Exact and Large Sample ML Techniques for Parameter
% Estimation and Detection in Array Processing.
%
% SEE ALSO
% wsf_detection.m, glrt_sml.m, glrt_experiment, apply_detection

% Copyright (C) 2018 Brian Emery
%
% Version 4-Jul-2018

% TO DO
% - code here will probably be called as part of a detection/estimation
% scheme that loops over n, and stops when the LR is small enough. It seems
% that 'small enough' is the criteria.


% OVSN93 DEFINITIONS
% m sensors
% d signals
% N snapshots (we use K)


% check for test case
if strcmp('--t',A), test_case, return, end



% # elements x # signals
[m,d] = size(A);

% Get identity matrix of correct size
I = eye(m);



% GET THE MODEL COV
% Here I'm useing O93 to construct the model cov based on the given
% solutions. This is how it's done in TF09, A.1.10, pg 48. Note that the
% input solutions could come from MUSIC or any MLE method (I suppose?).
% See also, the TF09 method --constructs cov matricies to compare without using
% the data cov, and then uses a threshold to determine which is more
% valid.--
% 

% Get the orthogonal projector onto the null space of S (4.41 of OVSN93)
Pp = I - A*pinv(A);

% (re) compute the Vdml 
Vdml = trace( Pp*R );

% Compute noise power (4.49, note that this differs from 4.38 though which
% has (m-d) in denominator
po = 1./(m) * Vdml;

% Compute the signal power (non-dB, units must be volts^2?)
S = signal_power(A,R,I,po);

% % ... S should be very similar to this:
% s = D.s;
% (1/K)* s*s' % ... check


% Compute the model covariance
Rm = A*S*A' + po*I;

% Rm and R should be pretty close if m is correct
%
% 
% also, this is an est of SNR I think:
% e = D.e;
% 10*log10( trace((1/K)* s*s') / trace( (1/K)* e*e' ) )



% Comptue the Likelihood Ratio for the given parameters 

% Note that Ch 7 of TF09 has arg = R*inv(Rm), which produces similar
% results in the tests at least (my Rm is the R_model)
arg =  R/Rm; % R*inv(Rm);
% arg = Rm\R; %inv(Rm)*R; 

% 2006 book chapter, equation 2.27
% ... also TF09 eqn 7.15. Note that when Rm=R, LR = 1
LR = ( det(arg)/ (trace(arg)/m)^m )^K;  %keyboard

% output real part only - imag part very small
LR = real(LR);

% .. or part of eqn 2.6
gm = (exp(1)/K)^(m*K) * ( det(arg)^K ) * exp(-trace(arg)) ;


% Output the simple LR from Ottersten 1993 (4.140 and 4.141)
%
% % they have Rm - R, but this makes sense i think? 
%
% Aug 2019:
% ... this is what 1998 swindlehurst et al have (eqn 81), with
% LR2 > threshold meaning ~model signal is correct, and LR < threshold
% meaning the model not inline with the signal (... or something)
LR2 = 2*K* ( log10(norm(R)) - log10(norm(Rm)) ); 

% TF09 eqn 7.11
% LR = ( (det(arg)*exp(m) )/exp( trace(arg) ) ).^K
 


% 
% % OVSN93 pg 50, step 5
% LR_5 = 2*N * ( Vdml - det(R)); %log(det(R)) ); ???
% %
% % ... or maybe this ...
% LR_5b = 2*N * ( log( det(S*P*S' + po*I) ) - log(det(R)) );
% 
% % see also the book with different weights etc
% 
% keyboard



end

 

function S = signal_power(A,R,I,po)
% SIGNAL POWER
%
% INPUT
% S(theta) at the DOAs
% 
% equations 8 and 9
% 
% See OVSN93 equation 4.126, and text therein which suggests using "the
% corresponding SML estimate (4.37) [for S]
% 
% Based on Schmidt 1986 this should generate the signal covariance matrix

% % *** 
% 
% % From the 2006 book chapter
% 
% % beta
% % B = inv(S'*S);
% 
% % noise power (mean of noise eigen values)?
% % po = 1/(M-m) * sum()
% 
% % P = (B*S') * [ R - po*eye(M)] * (S*B);
% 
% % the above, per matlab suggestion:
% P = ((S'*S)\S') * ( R - po*eye(M)) * (S/(S'*S));
%
% ... which looks like the SML estimate below ...
%
% % *** 


% 4.37 ... also similar to eqn 7 in Schmidt 1986 

Adag = pinv(A);

S = Adag * (R - po*I) * Adag';


end

function test_case
% TEST CASE
% 
% test case directory: /m_files/test_data/

% ok, trying this, just generate pdf's:
generate_pdfs





% % TF09 CH 7 CASE
% % "In Abramovich et al. (2007b), the absolute maximum LR(R0) 
% %  values associ- ated with a properly ordered model for a 12-element 
% %  antenna with 24 samples (i.e., a circumstance where the sample 
% %  covariance matrix is reasonably repre- sentative of the true 
% %  covariance matrix R0; Reed et al., 1974) were shown to average 
% %  around 0.03 and never exceeded 0.1, compared to the LR value of 
% %  unity associated with the unstructured MLE estimate (R? M ), 
% %  which has an LR of unity.
% th = -85:0.1:85;
% M = 12;          % # elements in array
% K = 24;          % # snapshots
% SNR = 20;
% 
% [A,~] = make_ula_pattern(th,M);
% 
% [~,ix,~] = intersect(th,0); 
% 
% 
% [R,~] = radar_simulation_basic( A(:,ix) ,K, SNR);
% 
% LR = glrt(A(:,ix),R,K); % = 6.5383e-46
% 
% keyboard


% BASIC VALIDATION TEST
% high SNR, high K, ... GLRT Assumes K > M (snapshots >
% antennas)
%
% Aug 2019: testing the Ottersten et al 1993 formulation

% ... of which this is similar, they correlate emitters though, and etc
th = -85:0.1:85;
M = 6; %8;          % # elements in array
K = 10;          % # snapshots
SNR = 6;

% Make the pattern (xu seems to work best?)
[A,~] = make_ula_pattern(th,M); %,'xu'); %,'moses'); or 'xu', all relative to perp

% [A,th] = make_ura_minus_one; % 78-320

[~,ix,~] = intersect(th,[0 10]); %[90 120 150]); % 1 2 3 and 4 ... 

figure

for i = 1:10
    
    [R,~] = radar_simulation_basic( A(:,ix) ,K, SNR);
    
    
    % keyboard
    %
    % % get ml estimates for various models ...
    % ... seems that the result is pretty consistent, around 1000 ...
    % ... with some probability of detection based on the snr ...
    % ... need to figure out how to make the empirical pdf ...
    
    idx = mle_ap(A,R,1);
    
    [LR(1),~,~,LR2(1)] = glrt(A(:,idx),R,K);
    
    % compute the chi^2 stat, c.f. https://www.mathworks.com/help/stats/chi2inv.html
    chi2(1) = chi2inv(0.05, dof(M,length(idx)) ); % 95% is arbitrary
    
        
    idx = mle_ap(A,R,2);
    
    [LR(2),~,~,LR2(2)] = glrt(A(:,idx),R,K);
    
    chi2(2) = chi2inv(0.05, dof(M,length(idx)) );
    
    
    
    idx = mle_ap(A,R,3);
    
    [LR(3),~,~,LR2(3)] = glrt(A(:,idx),R,K);
    
    chi2(3) = chi2inv(0.05, dof(M,length(idx)) );
    
    
    
    idx = mle_ap(A,R,4);
    
    [LR(4),~,~,LR2(4)] = glrt(A(:,idx),R,K);
    
    chi2(4) = chi2inv(0.05, dof(M,length(idx)) );
    
    keyboard
    
    % ... this below is older ... 
%     % hmm
%     plot(1:4, real(-2*(log(LR))) ) % not sure about 2*abs(log(LR))
%     
%     
%     hold on
%     
%     axis([1 4 900 1200])
    
    
    
end


end

function n = dof(m,d)
% DOF - compute degrees of freedome
% Ottersten 1993, pg 50

p = 1; % we're only looking for one parameter (direction)

n = (m.^2) - (d.^2) - (p.*d) - 1;




end


function previous_to_previous



% ----  PREVIOUS CODE ---

% set number of parameters (see below)
p = 1; 


% LOOP OVER POSSIBLE NUMBERS OF SIGNALS

% for testing, I want to look at the ratio, and the gamma computed
[gamma,LR ] = deal(NaN(1,m));


% set initial number of signal sources. They start with zero, which wtf?
for d = 1:m-1
    
    
    % Compute Gamma
    %
    % Find a value that exceeds 95% of the samples from a chi-square distribution with 10 degrees of freedom.
    %
    % x = chi2inv(0.95,10)
    % x =
    %   18.3070
    %
    % You would observe values greater than 18.3 only 5% of the time by chance.
    % gamma = 1; % NOTE, this is a function of the degrees of freedom, which is
    % dependent on the number of emitters
    
    % Compute initial estimate of the degrees of freedom
    % m^2 - d^2 - pd -1 in the paper, with d signals, m sensors, and p
    % parameters. For us, p = 1, since we're only looking for the DOA (see pg 9
    % of the reference).
    dof = m.^2 - d.^2 - p.*d - 1;
    
    % Compute the threshold based on the 95% value of the X^2 (chi squared)
    % distribution for this value of the degrees of freedom:
    gamma(d) = chi2inv(0.95,dof);
    
    
       
    % Use SML to find the do under the hypothesis that the number of
    % sources is d
    [ix,Vsml] = sml_ap(A,R,d);

    
    % --Evaluate the hypothesis test (goodness of fit)--
    % compute the ratio .. .might need real and imag parts of Vsml
    LR(d) = 2 * K * ( Vsml - real(log(det(R))) ); % maybe need just real part
    
    
    
end


keyboard



% prevous to delete


% initialize the loop
Vsml = 100;

d = 1;

%
% supposed to have the 2K here ...
% while ( 2*K*(Vsml - log(real(det(R)))) > gamma && d < N)
while ( (Vsml - log(real(det(R)))) > gamma && d < N)

    disp(['d = ' num2str(d) '... '])
    
    % increment
    d = d + 1;
    
    
end








end


function generate_pdfs
% % try to make the LR pdfs, for N = 1, 2, etc
% 
%  BASIC VALIDATION TEST
% high SNR, high K, ... GLRT Assumes K > M (snapshots >
% antennas)
%
% Aug 2019: testing the Ottersten et al 1993 formulation

N1 = run_sim(0);

N2 = run_sim([0 10]);

N3 = run_sim([-10 0 10]);

keyboard


end

function N = run_sim(src)
% RUN SIM - with input source bearings
%
% ... of which this is similar, they correlate emitters though, and etc
th = -85:0.1:85;
M = 6; %8;          % # elements in array
K = 10;          % # snapshots
SNR = 25*rand(100,1);

% Make the pattern (xu seems to work best?)
[A,~] = make_ula_pattern(th,M); %,'xu'); %,'moses'); or 'xu', all relative to perp

% [A,th] = make_ura_minus_one; % 78-320

[~,ix,~] = intersect(th,src); %[90 120 150]); % 1 2 3 and 4 ... 

[LR,LR2] = deal(nan(size(SNR)));

N.LR = LR;
N.LR2 = LR2;

for i = 1:numel(SNR)
    
    [R,~] = radar_simulation_basic( A(:,ix) ,K, SNR(i));
    
    
    % keyboard
    %
    % % get ml estimates for various models ...
    % ... seems that the result is pretty consistent, around 1000 ...
    % ... with some probability of detection based on the snr ...
    % ... need to figure out how to make the empirical pdf ...
    
    idx = mle_ap(A,R,length(src));
    
    [N.LR(i),~,~,N.LR2(i)] = glrt(A(:,idx),R,K);
       
    
end



end
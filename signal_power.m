function S = signal_power(A,R) %I,po)
% SIGNAL POWER - from glrt.m ... 
%
% INPUTS
% A    - Array matrix at theta(s) given by the DOA solution,
%        # elements x # signals.
% R    - Data Covariance matrix
%
% OUTPUTS
% S(theta) at the DOAs
% 
% equations 8 and 9
% 
% See OVSN93 equation 4.126, and text therein which suggests using "the
% corresponding SML estimate (4.37) [for S]
% 
% Based on Schmidt 1986 this should generate the signal covariance matrix
%

% TO DO 
% - See also 2004 abramovich which as Ottersten 1993 based 'traditional
%   method' ... look for test data 
% - Use alt estimates of noise? 
% - output diagnal elements only? 
%

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

% check for test case
if strcmp('--t',A), test_case;  end


% # elements x # signals
[m,d] = size(A);

% Get identity matrix of correct size
I = eye(m);



% Get the orthogonal projector onto the null space of S (4.41 of OVSN93)
Pp = I - A*pinv(A);

% (re) compute the Vdml 
Vdml = trace( Pp*R );

% Compute noise power (4.49, note that this differs from 4.38 though which
% has (m-d) in denominator
po = 1./(m) * Vdml;


% 4.37 ... also similar to eqn 7 in Schmidt 1986 

Adag = pinv(A);

S = Adag * (R - po*I) * Adag';

% output just the powers, not the whole matrix?
S = diag(S);

% keyboard
% 
% % Schmit 1986 eqn 7 added in reprint, gets same as above
% So = po*I; 
% ASAi = inv(A'*inv(So)*A);
% S = ASAi*A'*inv(So)*(R - po*I)*inv(So)*A*ASAi; %  same!

end


function test_case
%
% TEST CASE
%

% TEST UNDERSTANDING OF SIGNAL POWER
% also, make sure I can model different powers, and then resolve them with
% both MUSIC and MLE

% Setup
th = -85:0.1:85;
M = 12;          % # elements in array (not actually used below) 
K = 100;          % # snapshots K > M ... really big for this test



% Make the pattern (xu seems to work best?)
[A,~] = make_ula_pattern(th,M,'xu'); %,'moses'); or 'xu', all relative to perp

% Source indecies
[~,ix,~] = intersect(th,[0 45]);

Am = A(:,ix);

d = length(ix); % number of sources



% BASED ON ULADATA.M ... and from radar_simulation_basic.m
% from http://www2.ece.ohio-state.edu/~randy/SAtext/
%
% See 6.2.21
SNR = 30; % dB

% get noise power (power, so we use the /10 and not /20 in conversion)
% ... assumes max signal is 1 and gets noise relative to this
sig2 = 10^(-SNR/10); % the noise variance, or std squared

% get signal correlation matrix 
P = eye(d);


% generate the source signals. 
% SAS-book has power as the expected value of (det(s))^2  (6.3.21) ... but
% this implies that s is square so ...
% (Sept 2019 notes: this has randn producing a zero mean normally 
% distributed random numbers with a standard deviation of 1 ... for two
% sources (or more) these will be uncorrelated if P is an identity matrix)
s = (sqrtm(P)')*randn(d,K);     

% 1989 Kelly suggests that this can be done, ie that I can model one signal
% as being lower power than the other (and still uncorrelated)
% ... also The power is the MSE, and since these are zero mean, bias = 0, 
% and thus MSE = var 
%
s(1,:) = randn(1,K); % std of 1
s(2,:) = sqrt(0.05)*randn(1,K); % var of 0.5



% See SAS-book 6.4.2: signal covariance matrix
Ps =  (s*s')./K; % 




% generate the noise component
e = sqrt(sig2/2) * (randn(M,K)+1i*randn(M,K));

% generate the ULA data (time series) 
Y = Am*s + e;


% generate the data covariance matrix
R =  (Y*Y')./K;



% COMPUTE DOAs
ixml = mle_ap(A,R,2);

[~,ixmu,~,~] = music(A,R,th,2);

ixmu = ixmu{2};


% These get the correct answer!
S = signal_power(A(:,ix),R); % ixmu, ixml


% % recover the snr in db?
% 10*log10(S(1,1)/sig2)
% 10*log10(S(2,2)/sig2)

% recover the snr in db?
10*log10(S(1)/sig2)
10*log10(S(2)/sig2)



keyboard



% --------------

% Basic Test that suggests this is working and is able to resolve signals
% with different powers

%   205   330
[C,APM] = test_data_tonys;

A = get_array_matrix(APM);

ix = mle_ap(A,C,2);

S = signal_power(A(:,ix),C)




end


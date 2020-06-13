function [Pvalues,testResult] = run_param_test(eigValues,eigVectors,A,dualIdx,P)
% RUN MUSIC PARAMETER TESTS.M - determine 1 or 2 bearing solution
% [Pvalues,testResult] = run_param_test(eigValues,eigVectors,A,dualIdx,P)
% 
% Apply Codar's MUSIC parameters to determine if the 2bearing or 1bearing
% hypothesis is more valid. Since the test is run on the dual brg solution,
% testResult is * TRUE IF DUAL *, false if single.
%
% INPUTS 
% eigValues,eigVectors - of the covariance matrix
% A                    - antenna pattern data matrix ( m antennas x nbearings )
% dualIdx              - two indicies of the APM bearings (dual brg soln)
% P                    - the MUSIC parameters, eg [20 10 3] or [40 20 2]
%
% OUTPUTS
% Pvalues              - Computed paramter values
% testResults          - boolean, true if DUAL BEARING
%
% References:
% SeaSonde_Radial_Processing.pdf
% Radar Angle Determination with MUSIC Direction Finding, Nov 1999, U.S. Patent 5,990,834
%
% NOTES
% % plot the parameters in 3d space with the cube of dual == true
% plot3(S.Params(:,1),S.Params(:,2),S.Params(:,3),'.'), hold on
% grid
% xlabel('P1'),ylabel('P2'),zlabel('P3')
% axis([0 100 0 100 0 20])
% h = patch([0 40 40 0 0],[0 0 20 20 0],[2 2 2 2 2],'c');
% h = patch([40 40 40 40 40],[0 20 20 0 0],[2 2 100 100 2 ],'c');
% h = patch([0 40 40 0 0],[20 20 20 20 20],[2 2 100 100 2 ],'c');
% 
%
% See also:  apply_test_result.m, detection_codar_post_proc.m 

% Brian Emery 3 Giugno 2008
% Notes:
% Run the check vs tony's paper to get the music param values. the numbers
% I get are slightly different from his. I'm not sure why. This mfile
% should be pretty solid. The method of combining the real and imag parts
% of the APM is the only possible suspect. (29Jul08)
%
% 29 Jan 2015
% made this into a function, added test case
%
% 16 Jun 2015
% Coded a solution to the problem of detecting two peaks in the music
% spectrum, where only one peak is detected, and dualIdx is numel = 1.




% check for test case
if strcmp('--t',eigValues), test_case; return, end

% check for null case
if any(isnan(dualIdx)) || length(dualIdx) <= 1
    
    Pvalues = NaN(1,3);
    testResult = false;
    
    return
end


% INPUTS

% Set music parameters? 
% see [[music setting notes]]. This allows more dual solutions
if nargin < 5, P = [40 20 2]; end

% input must be vector of eigenvalues
if ~isvector(eigValues), eigValues = diag(eigValues); end

% trivial check smallest to largest
[eigValues,ie] = sort(eigValues); eigVectors = eigVectors(:,ie); 

% initialize the actual values for output
Pvalues = NaN(1,3); 

% Get the APM manifold ... code below wants it nbearings x m antennas
A = A.';



% RUN TESTS
% TWO bearing solution if:

% 1) the ratio of the largest eigenvalue to the 2nd largest is less
% than P1 
Pvalues(1) = eigValues(end)./eigValues(end-1); 



% 2) "the ratio of the largest two signal powers to the smallest [of the
% two signal powers] is be less than P2. This based on COS's Patent
% Version.

G  = conj(A(dualIdx,:))*eigVectors(:,[3 2]);
Gt = eigVectors(:,[3 2])'*A(dualIdx,:)';

% S  = inv(Gt)*[eigValues(end) 0; 0 eigValues(end-1) ]*inv(G);
S  = Gt\[eigValues(end) 0; 0 eigValues(end-1) ]/(G);

% Get Signal Power Ratio:
sigPowers = sort(diag(S)); 

Pvalues(2) = real(sigPowers(end))./real(sigPowers(end-1));



% 3) for the signal matrix, the ratio of the product of the diagonal
% elements to the product of the off diagonal elements must be greater
% than P3.
Pvalues(3) = real(prod(diag(S)) ./ prod(S([2 3])) );

% if isnan(Pvalues(3)), keyboard, end

% Run test (TRUE IF DUAL BEARING)
if Pvalues(1) < P(1) && Pvalues(2) < P(2) && Pvalues(3) > P(3)
    testResult = true;
else
    testResult = false;
end



end


function test_case
% TEST CASE

% EXAMPLE FROM 
% Properties of HF RADAR Compact Antenna Arrays and Their Effect 
% on the MUSIC Algorithm by dePaolo and Terril

% Create the covariance matrix from de Paolo's example:
C = [ 0.2162          0.0303-0.0090i  0.3170-0.0063i; ...
      0.0303+0.0090i  0.0436         -0.0091+0.0213i; ...
      0.3170+0.0063i -0.0091-0.0213i  0.5416];

[V,D] = eig(C);


APM = make_ideal_pattern(225);
A = get_array_matrix(APM);

%[DOA,singleIdx,dualIdx] = getDOAs(A,V,D);
% instead of using getDOAs, use the fact that the known signal input
% bearing are 205 and 330 deg (from the reference)
dualIdx = [find(APM.BEAR == 205) find(APM.BEAR == 330)];


[Pvalues,testResult] = run_param_test(D,V,A,dualIdx,[20 10 3]);


% check with values given in figure 9
if isequal(round(Pvalues*100)/100,[11.2800    4.4300    2.7200])
    disp('test case ... ok')
else
    disp('test case ... NOT ok')
    keyboard
end

return

figure
subplot(2,1,1), plot(A.BEAR,10*log(DOA(:,1)))

subplot(2,1,2), plot(A.BEAR,10*log(DOA(:,2)),'-r')






end

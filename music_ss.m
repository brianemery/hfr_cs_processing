function [DOA,singleIdx,dualIdx,testResult,Pvalues] = music_ss(APM,C)
% MUSIC -- method for direction of arrival in SeaSonde processing
%
% Does the eigen decomposition and DOA function, then runs the 
% (single vs dual) bearing test
%
% INPUTS
% C    - covariance matrix
% APM  - APM structure
%
% OUTPUTS
% DOA        - The DOA 'spectrum' for both single and dual, vs bearing
% singleIdx  - Index of APM for single angle soln
% dualIdx    - Indecies of APM for dual angle soln, might be NaN
% testResult - testResult is true if dual, false if single. 
% Pvalues    - Results of parameter test (see run_param_test.m)
% 
% Was music.m, which is now being generalized. This here for bw
% compatibility

% Copyright (C) 2015 Brian Emery

% TO DO
% check that units are NOT dbm
% more compatibility/swappability with mle.m, mle_ap.m, etc
%
% Add dePaolo's metrics?
%
% FUTURE:
% - less seasonde-centric
% - use array matrix as input, with bearing (theta)
% - input the number of signals to look for, and output solutions for
% - decide on a standardized output format for the pre-radial data


% GET EIGEN VECTORS AND EIGEN VALUES
% note: "[V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
% full matrix V whose columns are the corresponding eigenvectors so
% that X*V = V*D.", ie C*eigVectors=eigVectors*eigValues
[eigVectors,eigValues] = eig(C);


% Right now this is handled in getDOAs.m, but could do this instead. MATLAB
% appears to sort smallest to largest, but it's not a documented feature
% 
% ...] = eig(...,'vector') returns eigenvalues in a column vector 
%     instead of a diagonal matrix.
% 
% Then ...
% [eigValues,ix] = sort(eigValues);
% eigVectors = eigVectors(:,ix);



% Compute DOA's
[DOA,singleIdx,dualIdx] = getDOAs(APM,eigVectors,eigValues);



% ADD LABEL ABOUT HYPOTHESIS TEST RESULT
% testResult is true if dual, false if single.
[Pvalues,testResult] = run_param_test(eigValues,eigVectors,APM,dualIdx);



end

function [DOA,singleIdx,dualIdx] = getDOAs(APM,eigVectors,eigValues)
% GETDOAS.M - Get Directions of arrival using MUSIC
% [DOA,singleIdx,dualIdx]=getDOAs(APM,eigVectors,eigValues)
%
% Compute the possible directions of arrival assuming single bearing
% and dual bearing solutions. Project the Antenna manifold on to the
% noise subspaces for each Bearing. The bearing with the smallest
% projection is the most orthogonal and is likely where the signal
% came from. (copy here is from the simulateDivergenceProblem.m
% subfunction)
%
% INPUTS
% Antenna Pattern (APM structure)
% eigenVectors from the signal covariance matrix, formed from a single
% cross spectra freq bin, in a single range cell (eigVectors should be smallest to largest)
%
% OUTPUTS
% the DOA function (inverse dist squared as function of bearing), for
% single bearing solutions (column 1) and dual bearing solutions (col 2)
%
% DOA typically converted to dBm using 10*log10(DOA)
%
% Notes:
% Ideal A(theta) (APM Amplitude) has dimensions 1x3, En is noise eigen
% vector(s) with dimensions (3x1 or 3x2).
% This outputs potential 1 and 2 bearing solutions.
%
% Generally: DOA=1/(A'*En*En'A);
%
% Using the outputs:
% brgs1=APM.BEAR(singleIdx);
% brgs2=APM.BEAR(dualIdx);
%
%
% Schmidt't paper calls this the inverse of the euclidean distance squared.
% Since A and En are likely complex, I break these into each component (real
% and imag parts are seperate components) then compute dist.
%
% REFERENCES
% Based on De Paolo and Terrill Scripps report, along with:
% Barrick and Lipa 1983 (IEEE Journal of Ocean Engineering)
% Radar Angle Determination with MUSIC Direction Finding, Nov 1999, U.S. Patent 5,990,834


% Copyright (C) 2008-2010 Brian M. Emery

% Brian Emery Mayo 2008
%  3 Nov 2011 
%   - added test case and did more checking vs Tony's paper.
%   - discovered the difference between ' and .' and implemented it (arg!)
%
% TO DO
% better method for determining 2nd peak for dual bearing solutions (need
% simulations for this? or ship data?)

if strcmp('--t',APM), test_case, end

% CHECK INPUTS
% eigVectors should be smallest to largest
% Since this is a new feature, throw a warning if 2 inputs used
if nargin < 3,
    warning(id,'getDOAs.m now requires 3 inputs')
else

    % Sort Eigenvalues and EigenVectors smallest to largest. Note
    % that this is the opposite of the way codar does it.
    %
    % diag outputs a column vector in this case
    if ismatrix(eigValues), eigValues = diag(eigValues); end
    
    % sort smallest to largest which appears to be the matlab default
    [~,edx] = sort(eigValues);
    eigVectors = eigVectors(:,edx); clear edx

    % Could do this instead ...
%        [...] = eig(...,'vector') returns eigenvalues in a column vector 
%     instead of a diagonal matrix.
% 
    
end
 

% Initialize the Output
DOA(1:length(APM.BEAR),1:2) = NaN;

% Noise space for 2 brg soln 
E2 = eigVectors(:,1);

% Noise space for 1 brg soln
E1 = eigVectors(:,1:2);
                 
% Note that:
% sqrt(sum(E1(:,1).^2)) == 1 (Tony's section 4.2)


for idx=1:length(APM.BEAR)
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
        
    % 1 bearing solutions 
    % NOTE: A' is the complex conjugate transpose
    % sizes:      1x3 * 3xm * mx3 * 3x1
    DOA(idx,1)=1/(A'*E1*E1'*A); 
    
    % 2 bearing solutions (hmmm ...)
    DOA(idx,2)=1/(A'*E2*E2'*A); 
    
end


% Get indicies
[y,singleIdx] = max(DOA(:,1));

% [y,dualIdx]=sort(DOA(:,2)); 

% TRY PEAKFINDER.M
% sel = (max(x0)-min(x0))/4)
% some of the settings may need tweaking:
% 
% Do this in log space
x = 10*log10(real(DOA(:,2)));
[dualIdx,y] = peakfinder(x, (max(x)-min(x))/100, [], 1,1);

if length(dualIdx) > 2
    
    % sometimes the peak picking outputs both ends and the peak
    [y,ii] = sort(y);
    y = y(2:end);
    dualIdx = dualIdx(ii(2:end));
    
end

% NOTE 2:
% outputting NaN index is sort of pointless. Just need later code (eg
% doa_on_range_cell.m) to be able to have 0, 1 or 2 indecies from this
% function
% 
% % NOTE: this needs to always output 2 indecies? or to propagate single brg
% % solution
% if length(dualIdx) <= 1,  
%     dualIdx = [dualIdx NaN(1,2-length(dualIdx))]; 
% end

% previous to peakfinder.m usage:
% I need a better method for determining the dual bearing indicies. The
% test case should get 205 deg and 330 deg, rather than the adjacent
% bearings on the main peak. Could define a search region excluding the
% main (first) peak, base on the resolution (fxn of SNR?)

% keyboard

% dualIdx=[dualIdx(end) dualIdx(end-1)];





end


function test_case
% TEST CASE
% 
% Makes a figure showing the DOA function vs bearing just like fig 9 in the
% De Paolo and Terrill Scripps report
%
% Lots of code from music.m

% (APM,DOA,singleIdx,dualIdx)

% Create the covariance matrix from de Paolo's example:
% NOTE that he made this from a simulation with currents input at 205 and
% 330 degrees which MUSIC gets kind of wrong
C=[ 0.2162          0.0303-0.0090i  0.3170-0.0063i; ...
    0.0303+0.0090i  0.0436         -0.0091+0.0213i; ...
    0.3170+0.0063i -0.0091-0.0213i  0.5416];

% Create the idealized pattern 
APM = make_ideal_pattern(225, 0:5:360);

% GET EIGEN VECTORS AND EIGEN VALUES
% note: "[V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
% full matrix V whose columns are the corresponding eigenvectors so
% that X*V = V*D.", ie C*eigVectors=eigVectors*eigValues
[eigVectors,eigValues]=eig(C);

% Compute DOA's
[DOA,singleIdx,dualIdx] = getDOAs(APM,eigVectors,eigValues);

% FIGURES
figure

subplot(211)
plot_doa(APM,DOA(:,1),singleIdx,'Single Bearing DOA function')
axis([0 360 -4 12])


subplot(212)
plot_doa(APM,DOA(:,2),dualIdx,'Dual Bearing DOA function')
axis([0 360 -5 30])


% DISPLAY CHECK NUMBERS
check_results(DOA(:,1),APM.BEAR,225)
check_results(DOA(:,2),APM.BEAR,205)
check_results(DOA(:,2),APM.BEAR,330)

keyboard

ix = find(APM.BEAR == 330)
DOA(ix,2)
10*log10(DOA(ix,2))
ix = find(APM.BEAR == 205)
10*log10(DOA(ix,2))

end

function plot_doa(APM,DOA,idx,titleStr)
%
% idx gives the location of the max of the DOA, single or dual angle
% solutions

% % Dual angles might have NaN index
% idx = idx(~isnan(idx));

plot(APM.BEAR,10*log10(DOA),'-b.')

hold on

plot(APM.BEAR(idx),10*log10(DOA(idx)),'g*')

xlabel('bearing (deg CWN)'),ylabel('10*log10(DOA)')

title(titleStr)

% add info to plot
for i = 1:length(idx)
text(APM.BEAR(idx(i))+10, 10*log10(DOA(idx(i))), ['(' num2str(APM.BEAR(idx(i))) ',' num2str(10*log10(DOA(idx(i)))) ')'])

end

end

function check_results(DOA,BEAR,brg)
% CHECK RESULTS
% 
% quickly check my resuls vs his

mine = num2str(real(10*log10(DOA(BEAR == brg))),4);

switch brg
    case 205
        his = '28.6';
        
    case 225
        his = '9.5';
        
    case 330
        his = '21.1';
end

disp(['DOA Metrics: his = ' his ' mine = ' mine ])


end
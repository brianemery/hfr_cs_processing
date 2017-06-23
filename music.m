function [DOA,idx,D,V] = music(A,C,th,n)
% MUSIC - Direction of Arrival from Multiple Signal Classification
%  [DOA,idx,D,V] = music(A,C,th,n)
%
% INPUTS
% A    - Array Manifold (aka array matrix at all theta (complex))
% C    - Covariance matrix
% th   - bearings associated with A
% n    - max number of signal sources (defaults to M-1)*
%
% *If a radar has M elements and theta has d bearings, A is Mxd
% 
%
% OUTPUTS
% DOA    - The DOA 'spectrum' for both single and dual, vs bearing
% idx    - Cell array with indecies of the n signal source solutions
%          eg idx{1} has the index of the single source soln, idx{2} is
%          dual, etc
% D      - Vector of eigenvalues
% V      - Matrix of eigenvectors
%
% 
% SEE ALSO
% From music_ss.m and others. Use run_parameter_test.m on outputs for 
% SeaSonde MUSIC parameters. See also music_error.m, music_error2.m

% Copyright (C) 2017 Brian Emery


% TO DO
% check that units are NOT dbm on C input
% more compatibility/swappability with mle.m, mle_ap.m, etc
%
% Incorp these in from R2017b:
% islocalmin and islocalmax Functions: Detect local minima and maxima in data
% minkand maxk Functions: Find the k smallest or largest elements in an arr


% Check for test case
if strcmp('--t',A), test_case, end

% Check for number of signals to search for
if nargin < 4, n = size(A,1)-1; end

% Get eigen-decomposition
[V,D] = eig(C); 

% sort smallest to largest which appears to be the matlab default but is
% not documented as such
[D,ix] = sort(diag(D));  
V = V(:,ix); 

% get number of elements
M = size(V,1);




% Initialize the Outputs
DOA(1:length(th),1:n) = NaN;

% loop over number of signals assumed
for j = 1:n
     
    % loop over bearings, computing the DOA for each of the n cases
    for i = 1:length(th)
        
        
        % Noise sub space for j brg soln
        % noise subspace is 1:(M-n), b/c the (M-n+1):M is the signal subspace
        En = V(:,1:(M-j)); 
        
        % Compute the DOA function (Schmidt 1986, equation 6)
        DOA(i,j) = 1/(A(:,i)'*(En*En')*A(:,i));
        
    end
    
end



% % USE FINDPEAKS
% % requires signal processing toolbox .. and testing too 
% %
% %    [...] = findpeaks(...,'NPeaks',NP) specifies the maximum number of peaks
% %     to be found. NP is an integer greater than zero. If not specified, all
% %     peaks are returned. Use this parameter in conjunction with setting the
% %     sort direction to 'descend' to return the NP largest peaks. (see
% %     'SortStr')
% % 
% 
% for i = 1:n
%     % just output the indicies
%     [~,idx{i}] = findpeaks(real(DOA(:,i)),'NPeaks',i,'SortStr','descend');
% end
% 
% keyboard

% Prealllocate
idx = cell(n,1);

% % Get indicies
[~,idx{1}] = max(DOA(:,1));

% Loop over remaining signals
for i = 2:n
        
    % ALT PEAKFINDER.M
    % sel = (max(x0)-min(x0))/4)
    % some of the settings may need tweaking:
    %
    % Do this in log space
    x = 10*log10(real(DOA(:,i)));
    
    % was dividing by 100
    [idx{i},y] = peakfinder(x, (max(x)-min(x))/200, [], 1,1);
    
    % check for and remove end points
         y =      y(~ismember(idx{i},[1 length(th)]));
    idx{i} = idx{i}(~ismember(idx{i},[1 length(th)]));
    
    % get largest n peaks
    if length(idx{i}) > i
               
        % plot(th,x), hold on
        % plot(th(idx{i}),x(idx{i}),'ro')
        % keyboard
        
        % get index of peaks to keep (k), then keep the i largest    
        [~,k] = sort(y);
        idx{i} = idx{i}(k(end+1-i:end));
                
    end
    
    
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
% 
% 
% [evar,ebias,mse] = deal('');


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

% Get the array matrix
A = get_array_matrix(APM); 

% Run the test
[DOA,idx] = music(A,C,APM.BEAR);




% FIGURES
figure

subplot(211)
plot_doa(APM,DOA(:,1),idx{1},'Single Bearing DOA function')
axis([0 360 -4 12])


subplot(212)
plot_doa(APM,DOA(:,2),idx{2},'Dual Bearing DOA function')
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
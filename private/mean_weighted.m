function abar = mean_weighted(a,wt)
% MEAN WEIGHTED - compute weighted mean
% 
% Computes the weighted mean of each row of a
% matrix after removing NaN's. 
% 
% INPUTS
% a, wt, values and weights. Weights can be a vector equal to the number of
% columns, otherwise must be a matrix the same size as 'a'
%

% TO DO
% need an example, needs to work on rows of matricies

% Copyright (C) 2009-2010 Brian M. Emery
% Jun 2009
% Jul 2012 - added test case, expanded for matricies

% INPUTS

% run test case?
if strcmp(a,'--t'), test_case, abar =[]; return, end

% check same number of columns
if size(a,2) ~= size(wt,2)
    error([mfilename ':TooFewWeightsColumns'])
end

% check number of rows
if size(a,1) ~= size(wt,1)
    
    % expand wt if this is true
    wt = ones(size(a,1),1)*wt;

end


% COMPUTE WEIGHTED MEAN

% set NaN's to zero in the data matrix and weights so that they do not
% contribute
i=isnan(a+wt);

a(i) = 0; wt(i) = 0;

% compute weighted mean
abar = dot(a,wt,2)./sum(wt,2);

% if we've inserted an entire row of NaN, set this back to NaN
j = find(sum(i,2) == size(a,2));

if ~isempty(j)
    abar(j) = NaN;
end



end

function test_case
% TEST CASE
% using data from wikipedia:


% VECTOR TEST
am = [62, 67, 71, 74, 76, 77, 78, 79, 79, 80, 80, 81,81, 82, 83, 84, 86, 89, 93, 98];

pm = [81, 82, 83, 84, 85, 86, 87, 87, 88, 88, 89, 89, 89, 90, 90, 90, 90, 91, 91, 91, 92, 92, 93, 93, 94, 95, 96, 97, 98, 99];

% straight mean
mn = mean_noNaN([am pm]);

% check by weight by n
am_ = mean_noNaN(am);
pm_ = mean_noNaN(pm);

% weight by n
wt = mean_weighted([am_ pm_],[length(am) length(pm)]);

% test
run_test(mn,wt)



% MATRIX TEST
dat = [am_ pm_; am_ pm_; am_ pm_;];
wt = mean_weighted(dat,[length(am) length(pm)]);

% test
run_test([mn mn mn]',wt)


% TEST WITH NANS
% want to make sure I'm not biasing the result, so the 2nd row should get 
% 86 (as before), different than having a zero 
dat(:,3) = [0 NaN 10]';

% add an entire row of NaN's also
dat(4,:) = NaN;

% define weights
wts = [length(am) length(pm) 30];

% compute
wt = mean_weighted(dat,wts);

% test vs pre-computed result
run_test([53.75 86 57.5 NaN]',wt)

end

function run_test(mn,wt)

if isequalwithequalnans(mn,wt)
    disp('testing mean_weighted.m ... ok')
else
    disp('testing mean_weighted.m ... not ok')
end

end



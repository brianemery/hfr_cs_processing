function R = apply_test_result(R)
% APPLY TEST RESULT - apply single/dual angle test results to DOA outputs
%
% Reform radial struct to contain only valid data (combine single and dual
% bearing results)

% Copyright (C) 2016 Brian Emery
% version 4 Apr 2016

if strcmp(R,'--t'), test_case, return, end

% Handle multi-element structs
if numel(R) > 1, R = struct_recursion(@apply_test_result,R); return, end
    

% list the fields to drop colums from? These are formatted such that
% columns are single dual dual 
fn = {'Bear','RadVel','dBear','mus_err','Err','mus_err2', ...
          'RngIdx','emitters','emitters_snr'};

% get only those that are present
fn = intersect(fn,fieldnames(R));

% Get single bearing solutions
S = subsref_struct(R,~R.Dual,size(R.Bear,1),1);

% Get rid of dual bearing data columns
% need to apply this to specific fields
S = drop_extra_cols(S,1,fn); 



% Get dual solutions
D = subsref_struct(R,logical(R.Dual),size(R.Bear,1),1);


% Get rid of single bearing data columns
% old way: D = subsref_struct(D,2:3,size(D.Bear,2),2);
D = drop_extra_cols(D,2:3,fn);



% Reshape the dual data: "elements are taken columnwise" by reshape
r = size(D.Bear,1);

for i = 1:numel(fn)
    D.(fn{i}) = reshape(D.(fn{i}),r*2,1);
end


% Vertically tile these 
fn = {'Params','Dual','SNR','RangeBearHead'}; 

for i = 1:numel(fn)
    D.(fn{i})= repmat(D.(fn{i}),2,1);
end


% now merge it all together 
R = struct_cat(1,S,D);

% clean up
R.TimeStamp =  unique(R.TimeStamp);
R.SiteOrigin = R.SiteOrigin(1,:);

end


function S = drop_extra_cols(S,col,fn)
% Like subsref struct, but for specific fields

% fn = {'Bear','RadVel'};

for i = 1:numel(fn)
    S.(fn{i}) = S.(fn{i})(:,col);
end


end

function test_case
% Basic test, not multi-element

load /m_files/test_data/apply_test_result.mat R

O = apply_test_result(R);


 keyboard


end

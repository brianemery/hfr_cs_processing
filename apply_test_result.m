function M = apply_test_result(R)
% APPLY TEST RESULT - apply single/dual angle test results to DOA outputs
%
% Reform radial struct to contain only valid data (combine single and dual
% bearing results)
%
% Compatible with the expanded DOA structs (see doa_struct.m)

% Copyright (C) 2016 Brian Emery
% version 4 Apr 2016

if strcmp(R,'--t'), test_case, return, end

% Handle multi-element structs
if numel(R) > 1, R = struct_recursion(@apply_test_result,R); return, end
    


% list the fields to drop colums from? These are formatted such that
% columns are single dual dual 
fn = {'Bear','RadVel','dBear','mus_err','Err','mus_err2', ...
          'emitters','emitters_snr'};

% get only those that are present
fn = intersect(fn,fieldnames(R));

% Get single bearing solutions
S = subsref_struct(R,~R.Dual,size(R.RadVel,1),1);

% Get dual solutions
D = subsref_struct(R,logical(R.Dual),size(R.RadVel,1),1);




% Get rid of dual bearing data columns
% need to apply this to specific fields
S = drop_extra_cols(S,1,fn); 

% Get rid of single bearing data columns
% old way: D = subsref_struct(D,2:3,size(D.Bear,2),2);
D = drop_extra_cols(D,2:3,fn);



% RESHAPE AND TILE DUAL DATA

% Reshape the dual data into column vectors: "elements are taken columnwise" by reshape
D = reshape_dual(D,fn);

% Vertically tile these ... if present
fn = intersect({'Params','Dual','SNR','RangeBearHead', ...
                                'RngIdx','eigValues'},fieldnames(R)); 

for i = 1:numel(fn)
    D.(fn{i})= repmat(D.(fn{i}),2,1);
end




% now merge it all together 
M = struct_cat(1,S,D);

% clean up
M.TimeStamp =  unique(M.TimeStamp);
M.SiteOrigin = R.SiteOrigin(1,:);

% no need to double these
if isfield(S,'RunTime')
    M.RunTime = S.RunTime;
end


end

function D = reshape_dual(D,fn,r)
% Reshape the dual data into column vectors: 
% "elements are taken columnwise" by reshape

if nargin < 3
    r = 2*size(D.RadVel,1);
end

for i = 1:numel(fn)
    
    if isstruct(D.(fn{i}))
        
        sfn = fieldnames(D.(fn{i}));
        
        D.(fn{i}) = reshape_dual(D.(fn{i}),sfn,r);
        
    else
        
        D.(fn{i}) = reshape(D.(fn{i}),r,1);
    end
end


end

function S = drop_extra_cols(S,col,fn)
% Like subsref struct, but for specific fields

% fn = {'Bear','RadVel'};

for i = 1:numel(fn)
    
    if ~isstruct(S.(fn{i}))
        
        S.(fn{i}) = S.(fn{i})(:,col);
        
    else
        
        % the field name refers to a struct, so get it's fields and apply
        % the indexing to these using recursion
        sfn = fieldnames(S.(fn{i}));
        
        S.(fn{i}) = drop_extra_cols(S.(fn{i}),col,sfn);
        
    end
end


end

function test_case
% Basic test, not multi-element

load /m_files/test_data/apply_test_result.mat R

S = apply_test_result(R);


 keyboard


end

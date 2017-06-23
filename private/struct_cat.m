function S = struct_cat(dim,varargin)
% STRUCT CAT - concatenate common fields in a struct using cat
% S = struct_cat(dim,S1,S2, ...)
% 
% Performs row or column concatenation of vector and matrix fields 
% (with the same name) from different structures. Unique is run on 
% strings and cells (these are assumed to contain meta data).
%
% INPUT
% similar, separate strucutures (common fields), or a single
% multi-element structure.
%
% EXAMPLE
% initialize a test struct
% S.a =  1:10;
% S.b = (1:10)';
% S.c = magic(5);
% S.str = 'strings';
% S.cll = {'cell strings'};
% 
% % create multi element struct
% M(1) = S; M(2) = S;
%
% % Row concatenate the fields:
% C = struct_cat(1,M);
%
% % Column concatenate seperate structures:
% C = struct_cat(2,S,S,S);
%
% see also field_concat.m, tuv_cat.m
%
% NOTE this may cause lots of side effects for the char fields 
% when dealing with radial structs


% Copyright (C) 2012 Brian Emery

% TO DO
% preallocated space?

% optionally run test
if strcmp('--t',dim), test_case, return, end

% get contents of varargin into multi-element struct
M = smart_struct_cat(varargin{:});

% M = [varargin{:}];

% get field names
fn = fieldnames(M);

% turn of warning related to use of 'rows' in unique below
warning('OFF','MATLAB:CELL:UNIQUE:RowsFlagIgnored')


% Run the concatenation
for i = 1:numel(fn)
    
    % row concat only for string and cell arrays
    if ischar(M(1).(fn{i}))
        
        % need to get rid of empties
        c = {M.(fn{i})}; c = c(~isempty(c));
        
        % now find the uniques
        S.(fn{i}) = unique(char(c),'rows');
        
    elseif iscell(M(1).(fn{i}))
        
        % force cells to row concatenate
        S.(fn{i}) = cat(1,M.(fn{i}));
        
    elseif isstruct(M(1).(fn{i})) & ~isempty( M(1).(fn{i}) )
        
        % handle possbily dissimilar sub structs with recursion
        S.(fn{i}) = struct_cat(dim,M.(fn{i}));
        
    else
        % arrays and matricies
        S.(fn{i}) = cat(dim,M.(fn{i}));
        
    end
    
end


end

function M = smart_struct_cat(varargin)
% SMART STRUCT CAT - allows for structs with different field names
% 
% code finds missing fields, add them to the
% other struct, then tries again
    

% Try the simple way (all fields named)
try 
    M = [varargin{:}];
    
catch
    
    M = varargin{1};
    
    for i = 2:numel(varargin)
       
        try 
            M(i) = varargin{i};
            
        catch %E
            
            % get dummy struct
            S = varargin{i};
            
            % find and add missing fields, both ways
            [M,S] = add_missing_fields(M,S);
            [S,M] = add_missing_fields(S,M);
            
            
            % now the concat should work. Allow for possibility of 
            % multi-element structure S
            M(i:i+(numel(S)-1)) = S;
            
            
            
        end
    end   
end
         
end

function [M,S] = add_missing_fields(M,S)
% ADD MISSING FIELDS

% find field in M that is not in S
fn = setdiff(fieldnames(M),fieldnames(S));

if ~isempty(fn)
    
    % add it to S
    for j = 1:numel(fn)
        [S.(fn{j})] = deal([]);
    end
    
    % SNEW = ORDERFIELDS(S1, S2) orders the fields in S1 so the new structure
    %  array SNEW has field names in the same order as those in S2. Sl and S2
    %  must have the same fields.
    S = orderfields(S,M);
    
end

end


function test_case 
% TEST CASE 

% handle dfferent fields? 

% init test struct
S.a =  1:10;
S.b = (1:10)';
S.c = magic(5);
S.d = (1:10)';
S.str = 'strings';
S.cll = {'cell strings'};
S.LAB = S;

% create multi element struct
M(1) = S; M(2) = S;

% make one empty
M(2).d =[];

% add extra field to one
S.newf = {'stuff'};

% Run Tests for inspection, output a 'report'
disp('S struct:')
S

disp('multiple input case (row concat)')
disp('S1 = struct_cat(1,S,S,S)')
S1 = struct_cat(1,S,S,S)

disp('multiple input case (col concat)')
disp('S2 = struct_cat(2,S,S)')
S2 = struct_cat(2,S,S)

disp('multi element input case, row concat')
disp('S3 = struct_cat(1,M)')
S3 = struct_cat(1,M)

disp('multi element input case, dissimlar fields')
disp('S3 = struct_cat(1,S,M)')
S4 = struct_cat(1,S,M)

keyboard


end
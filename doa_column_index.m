function col = doa_column_index(max_n)
% CREATE DOA COLUMN INDEX - map of number of emitters to the column index
%
% DOA structs have a matrix of (for example) the DOA solution for 1
% emitter, 2 emitter, etc solutions, where these are placed in columns.
% This creates the indexing that maps the number of emitters to the columns of
% the DOA matrix. (Note that the DOA stuct holds data prior to the 
% application of detection methods).
%
% EXAMPLE OUTPUT`
% col{:}
% ans =     1
% ans =     2     3
% ans =     4     5     6
% ans =     7     8     9    10
% ans =    11    12    13    14    15
% ans =    16    17    18    19    20    21
% ans =    22    23    24    25    26    27    28
%
% SEE ALSO
% doa_on_range_cell.m

% Copyright (C) 2020 Brian M. Emery

% ... worked out empirically ...

% preallocate
col = cell(1,max_n);

% create the indexing 
for i = 1:max_n
    
    % get the starting index
    i1 = sum(1:i) - i + 1;
    
    % create the column indexing
    col{i} = i1: (i1+i-1);
    
end

end


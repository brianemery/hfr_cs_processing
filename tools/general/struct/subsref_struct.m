function S = subsref_struct(S,idx,n,rc,fn)
% SUBSREF STRUCT - apply indexing to fields in a structure
% S = subsref_struct(S,i,n,rc,fn)
%
% Applies the indexing given in i, of the fields in S that have the 
% matching number of columns (n). Given a 4th input will apply the
% indexing to the rows (1). Defaults to columns (rc = 2). Cell arrays
% are now included, and input index array can be logical. 5th input can
% be a cell list of fields to skip over.
%
% Recurses into substructures looking for fields of the same size to apply
% the indexing to as well. 
%
% Example:
% % Create a structure with mixed fields, array and structure
% d.lon = ones(5,1)*[1:10]; 
% d.lat = ones(5,1)*[1:10]; 
% d.time = 1:10;
% d.strct = d; 
% 
% % Change the column indexing
% d = subsref_struct(d,2:5,size(d,2),2)
% 
% % Change the row indexing
% d = subsref_struct(d,1:2,size(d,1),1)
%
% SEE ALSO
% work_fields.m

% Copyright (C) 2010 Brian M. Emery
% Brian Emery 12 Jan 2010
% from subsrefTUV.m

% TO DO
% - maybe a special case when indexing applies to whole matrix ...?
% - should apply to cells also
% - recursion?
% - allow 3 inputs, last input is a size(X) and apply indexing to everyting
% size(X)
% - use structfun? might be a lot faster

if strcmp(S,'--t'), test_case, return, end

% don't skip any by default
if nargin < 5
    fn = '';
end

% Apply to columns by default
if nargin < 4 
    rc = 2; 
end

% Sort out indexing for each case
if rc == 2;
    rows = ':'; cols = idx;
elseif rc == 1;
    rows = idx; cols =':';
end

% % Recurse if multi-element
% if numel(S)>1, S = struct_recursion(@subsref_struct,S,idx,n,rc); return, end


% % < --- possible future case with absolute indexing
% elseif rc == 0 
%     [i,k] = ind2sub(size(),i); % <----


% Get field names to loop over
nm = fieldnames(S); 

% % experiment to make this faster, get only the field names we want to
% % operate on ... actualy makes it slower. The anon function seems to be the
% % problem
% sz = structfun(@(x) ( size(x,rc) ),S);
% 
% nm = nm(sz == n);

% detect char and struct arrays and skip them
s1 = structfun(@ischar,S);
s2 = structfun(@isstruct,S);

nm = nm(~s1 & ~s2);

% skim any named fields
nm = setdiff(nm,fn); 


for j = 1:numel(nm)
    
    %      keyboard
    %
    % Apply indexing to numeric or cell fields matching the correct number
    % of elements (rows or columns)
    if  size(S.(nm{j}),rc) == n  % && ( isnumeric(S.(nm{j})) || iscell(S.(nm{j})) || islogical(S.(nm{j})) )
        
        S.(nm{j}) = S.(nm{j})(rows,cols);
        
        
        %     % Recurse into sub-structures if they are numel == 1
        %     elseif isstruct(S.(nm{j})) && numel(S.(nm{j})) == 1
        %
        %         S.(nm{j}) = subsref_struct(S.(nm{j}),idx,n,rc);
        %
    end
end

% experitment: recurse here
if any(s2)
    nm = fieldnames(S); % list was shortened above
    nm = nm(s2);

    for i = 1:numel(nm)
        if numel(S.(nm{i})) == 1
            S.(nm{i}) = subsref_struct(S.(nm{i}),idx,n,rc,fn);
        end
    end
    
end
    
    

% % Optionally add documentation
% if isfield(S,'ProcessingSteps')
%     S.ProcessingSteps{end+1} = mfilename;
% end


end
%  ------------------------
function test_case
% Other tests needed? pass fail test ...?
% 

% % RECURSION TEST
% load /Data/testData/subsref_struct.mat
% 
% [~,t,j]=intersect(RDL.TimeStamp,round(DRFT.TimeStamp.*24)./24); clear j
% 
% R = subsref_struct(RDL,t,size(RDL.TimeStamp,2),2);
% keyboard


% ROBUST TEST
% Example FOR MATHWORKS :
% Create a structure with mixed fields, array and structure
d.lon = ones(7,1)*(1:10); 
d.lat = ones(7,1)*(1:10); 
d.time = 1:10;
d.strct = d; 
d.ProcessingSteps = {'blahblah'};
d.cellTst = cellstr(num2str([1:7]'))

disp(' Change the column indexing')
f = subsref_struct(d,2:5,10,2)

disp(' Change the column indexing')
f = subsref_struct(d,2:5,10)

disp(' Change the row indexing')
f = subsref_struct(d,[5 7],7,1)

keyboard

end

% ----------------------------------------------------------
% Additionally ...
function B = subsref_local(B,r,c,flds)
% SUBSREF LOCAL - local version of subscript reference
% Add this as a subfunction to clean up all the loops over fields. This
% applies the indexing in r, c to all fields listed in flds cell.
%
% at 181:
% % RDL = subsref_local(RDL,cat(2,sflds,flds),r,c);
%
% Note r or c can be ':' for all rows or columns.

if nargin < 4
    flds = fieldnames(B);
end

for i = 1:numel(flds)
    B.(flds{i}) = B.(flds{i})(r,c);
end

end


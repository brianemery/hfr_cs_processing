function field_check(D,flds)
% FIELD CHECK - check for presence of needed fields in input structures
% field_check(D,flds)
% 
% Stops mfile execution and returns an error if needed fields are not 
% found. Designed to be called by functions with structure inputs.

% Copyright (C) 2010 Brian Emery
% 2 Mar 10

% [st,i] = dbstack

% Check for critical fields
if any(~isfield(D,flds))
    error('Critical Structure Fields not found')    
end

end

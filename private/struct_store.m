function [APM,pl] = struct_store(APM,pl,S,keyfld)
% STRUCT STORE - index struct into pre-allocated struct with same fields
%
% [APM,pl] = struct_store(APM,pl,S,keyfld)
%
% Designed to be used in loops which create data stored in structures.
% Simplifies adding the data structure to a big, preallocated data
% structure.
%
% The structure to insert (eg S) can have a subset of the fields in APM.
% This function loops over the fields in S.
%
% INPUTS
% APM    - the pre-allocated structure
% pl     - the place marker (in this case the row index - 1, ie starts at 0)
% keyfld - (optional) fields matching this # rows are stored ('BEAR' default)
% S      - the struct to insert
%       
% ** Recurses into sub-structures with 'xbar' set as the key field **
%
% EXAMPLE:
% pl = 0;
%
% % Run the loop
% for i = 1:10
%     
%     % Create the struct data to insert
%     S = apm_struct(1);
%     
%     % Populate fields
%     for j = 1:numel(vars)
%         S.(vars{j}) = 2^i; 
%     end
% 
%     % Test substructs and multi-columns
%     S.RadVel.stdev = 2^i;
%     S.SNR = [2^i 2^i 2^i];
%     
%     
%     % Insert into the preallcated struct
%     [APM,pl] = struct_store(APM,pl,S); 
%     
% end
%
% % from the optional test:
% [APM,pl] = struct_store('--t');
%
% Also, to get rid of empties:
% i = find(~isnan(APM.A13R));
% APM = subsref_struct(APM,i,size(APM.A13R,1),1);


% Copyright (C) Brian M. Emery
% 15 November 2010

% SEE ALSO 
% subfunction to compute_proximity.m which uses index instead of place
% marker

% Optional test
if strcmp(APM,'--t'), test_case, [APM,pl] = deal([]); return, end

% use # rows of key field to control which fields are stored
if nargin < 4
    keyfld = 'BEAR';
end

% get # rows of key field. The try allows recursion below, catching only
% structs with the specified key field
try s = size(APM.(keyfld),1); catch return; end


% Get all the field names
flds=fieldnames(S);

% get size of field for inserting
sz = size(S.(keyfld),1);

for i=1:numel(flds)
 

    if s == size(APM.(flds{i}),1) && sz == size(S.(flds{i}),1)      
                
        % insert data
        APM.(flds{i})(pl+(1:sz),:) = S.(flds{i});
        
    elseif isstruct(S.(flds{i}))
        
        % recurse into structures. Here I'm requiring the substruct to key
        % off of the 'xbar' field
        [APM.(flds{i}),x] = struct_store(APM.(flds{i}),pl,S.(flds{i}),'xbar');
         
    end

end

% Advance the place marker
pl = pl + sz;


end
% ------------------------------------------------------------------------
function test_case

% Create preallocated struct to add to
APM = apm_struct(10);

% Init the place marker
pl = 0;

% vars
vars = {'BEAR','A13R','A13I','A23R','A23I','A13M','A13P','A23M','A23P'};

% Run the loop
for i = 1:10
    
    % Create the struct data to insert
    S = apm_struct(1);
    
    % Populate fields
    for j = 1:numel(vars)
        S.(vars{j}) = 2^i; 
    end

    % Test substructs and multi-columns
    S.RadVel.stdev = 2^i;
    S.SNR = [2^i 2^i 2^i];
    
    
    % Insert into the preallcated struct
    [APM,pl] = struct_store(APM,pl,S); 
    
end

keyboard

% APM.SNR
% APM.BEAR
% APM.RadVel.stdev

end
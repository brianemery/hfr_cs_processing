function M = struct_pack(fn,M)
% STRUCT PACK - move varibles into structure
% S = struct_pack(varibles, S)
%
% INPUT
% cell array with names of variable 
% (optionally) the struct to append to 
%
% EXAMPLE
% B = who;
% B = struct_pack(setdiff({B.name},'S'));
%
% - or -
% B = struct_pack(who);
%
% see also: struct_unpack, substruct_unpack

% Copyright (C) 2010 Brian M. Emery

if nargin < 2, M = []; 

else
    fn = setdiff(fn,inputname(2));
end 

for i = 1:numel(fn)
    try
        M.(fn{i}) = evalin('caller',fn{i});
    catch
        disp([mfilename ': ' fn{i} ' not packed: not found?'])
    end
end


end
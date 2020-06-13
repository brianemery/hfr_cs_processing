function struct_unpack(varargin)
% STRUCT UNPACK - move structure contents to workspace as variables
% struct_unpack(R)
% Allow for multiple structure inputs, eg struct_unpack(S,R,MET), etc.
%
% see also: struct_pack, substruct_unpack

% Copyright (C) 2010 Brian M. Emery

% NOTE
% ASSIGNIN Assign variable in workspace.
%     ASSIGNIN(WS,'name',V) assigns the variable 'name' in the
%     workspace WS the value V.  WS can be one of 'caller' or 'base'. 

for j = 1:numel(varargin)

    % get the structure
    R = varargin{j};

    % get the field names
    try
        fn = fieldnames(R);
        
        if ~isempty(fn)
            for i = 1:numel(fn)
                assignin('caller',fn{i},R.(fn{i}))
            end
        end
        
    catch
    end
end
end
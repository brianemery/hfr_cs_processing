function H = struct_recursion(fhdl,varargin)
% STRUCT RECURSION - apply mfile to each element of a structure
%
% Handles the recursion into structs (S) with numel(S) > 1, for mfiles that
% typically apply to the fields only.
%
% EXAMPLE
% see cs_volts2dbm.m, whos contents are applied to the fields, this line of
% code:
% 
% % OLD WAY - STILL USE THIS
% if numel(CS)>1, CS = struct_recursion(@cs_volts2dbm,CS); return, end
% 
% % NEW WAY (doesn't work ...)
% % Handles the if numel part internally:
% CS = struct_recursion(@cs_volts2dbm,CS);
% 
% makes sure that the mfile is applied to each element of CS

% TO DO
% needs work


S = varargin{1};

% recurse if multi-element MT input (how to sub-funct this?)
%if numel(S) > 1
    for i = 1:numel(S)
        H(i) = feval(fhdl,S(i)); %,varargin{2:end});
    end
% end

end
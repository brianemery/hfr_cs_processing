function save_ctf_formats(N,structStr)
% SAVE CTF FORMATS - merge new data with existing and save it
% 
% Called by ctf_formats.m and ctf_descriptions.m to append field 
% data to existing structures.
% 
% INPUT is the string field name to append to
%
% Quick and not particularly elegant

% Copyright (C) 2012 Brian Emery


% get ctf_formats.mat path
fname = which('ctf_formats.mat');

% load desc struct
S = load(fname,structStr);

% merge structs
M = [fieldnames(N) struct2cell(N); fieldnames(S.(structStr)) struct2cell(S.(structStr))];

% find unique fields, favoring the old over the new. 
% Using non-cell M, and unique(..,'rows') would be able to identify cases
% where keys have different descriptions
[~, r] = unique(M(:,1));
M = M(r,:); 

% set struct to save
eval([structStr ' = cell2struct(M(:,2), M(:,1), 1);'])


% write data (maybe just fprintf, as a .m?)
save(fname,structStr,'-append')


end


function extra_code
% EG ELIMINATE FIELDS

fieldnames(fmt)'

fn = {'SBG1'    'SBG2'    'SBG3'    'SIR1'    'SIR2'    'SIR3'    'SLD1'    'SLD2'    'SLD3' };

desc = rmfield(desc,fn);
fmt = rmfield(fmt,fn);
units = rmfield(units,fn);

end
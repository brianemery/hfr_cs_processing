function fname = fileparts_name(fullcsq)
% FILEPARTS NAME - get names only from list of full-path file names
% inputs are generally a cell with full path to csq names
%
% see also fileparts_cell.m

% TO DO 
% use cellfun!

% check ins, 
if ischar(fullcsq)
    fullcsq = cellstr(fullcsq);
end


fname = cell(numel(fullcsq),1);

for i = 1:numel(fullcsq)
    
    [~,fname{i},~] = fileparts(fullcsq{i});
    
end

% this shouldn't do uniques
% fname = unique(fname);

end

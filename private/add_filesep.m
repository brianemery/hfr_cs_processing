function wd = add_filesep(wd)
% ADD FILESEP - Detect trailing slash, insert if missing
% wd = add_filesep(wd)
% 
% char or cell input allowed

% TO DO
% rename this filesep_add
% make filesep_rm (see get_file_list)

% Copyright(C) 2011 Brian Emery

if strcmp('--t',wd), test_case, return, end

if ischar(wd)
    
    if ~strcmp(wd(end),filesep), wd(end+1) = filesep; end
    
elseif iscell(wd)
    
    for i = 1:numel(wd)
        
        if ~strcmp(wd{i}(end),filesep), wd{i}(end+1) = filesep; end
        
    end
    
end

end

function test_case

disp(['all should end in ' filesep])
add_filesep('wft')
add_filesep('wft/')
add_filesep({'wft/'})
add_filesep({'wft'})

keyboard

end
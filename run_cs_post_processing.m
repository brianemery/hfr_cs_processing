function run_cs_post_processing(wd,site)
% RUN CS POST PROCESSING - post process run_cs_processing.m radials
% 
% Concat structs, and run temporal averaging of the sub-hourly data
% Also does some clean up and parsing of fields to enable the temporal
% concatenation. 
%
% FROM
% gather_radial_data.m ... then the subfunction ... 

% Copyright(C) 2017 Brian Emery



% Process MUSIC derived radials
load_cat_save('Rmu',wd,site)


% Process MLE derived radials
load_cat_save('Rml',wd,site)



end

function load_cat_save(rfn,wd,site)
% LOAD CAT SAVE - load, concat and save data
%
% INPUT
% rfn - field name in the saved radial mat file (Currenlty either Rmu or
% Rml)
% wd - working directory (with RDLs subdirectory)


disp(['PROCESSING ' upper(rfn) ' DATA ...'])



% SETTINGS

% Fields to concat
fn = {'RadComp','SNR','Err','Dual','Params1','Params2','Params3','SNR3','SNR2'};


% GET RDL SHORTS
% wd = ['/home/emery/data/' site];
% 
flist = get_file_list([wd '/RDLs/'],'RDLm*');    % ****** % LIMIT <---------------
% flist = flist(1:240);

% get name of temporal cat'd radial struct
tfn = regexprep(rfn,'R','T');


% CONCAT DOA PROCESSED

for i = 1:numel(flist)
    
    S = load(flist{i});
    
    S.(rfn) = fix_multicol_fields(S.(rfn));
    
    R(i) = S.(rfn);
    
end

disp('... done loading')


eval([rfn ' = temporalConcatRadials_exact(R,fn);'])

disp('... done temporal concat ... doing hourly merge')



% Make into hourly data
eval([tfn ' = hourly_merge(' rfn ');'])

disp('... done merging concat')

save([wd '/radials_' site '_' tfn '_' datestr(now,'yyyymmdd') '.mat'],tfn,'-v7.3')
save([wd '/radials_' site '_' rfn '_' datestr(now,'yyyymmdd') '.mat'],rfn,'-v7.3')

disp('... done saving')



end

function R = fix_multicol_fields(R)
% FIX MULTICOL FIELDS - dole out fields with multi columns for cat in time


% Fix some fields prior to the temporal concat

R.Params1 = R.Params(:,1);
R.Params2 = R.Params(:,2);
R.Params3 = R.Params(:,3);

R.SNR3 = R.SNR(:,3);
R.SNR2 = R.SNR(:,2);
R.SNR = R.SNR(:,1);

R = rmfield(R,'Params');

R.FileName = {R.FileName};

end



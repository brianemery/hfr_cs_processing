function fn = cs_fieldnames(CS)
% CS FIELDNAMES - read antenna field names from Cross Spectra File
% fn = cs_fieldnames(CS)
% 
%     fn = {'antenna1Self'
%         'antenna2Self'
%         'antenna3Self'
%         'antenna12CrossSp'
%         'antenna13CrossSp'
%         'antenna23CrossSp'
%                           };
%
% SEE ALSO
% ... for generating them use this: cs_make_field_names.m


fn = fieldnames(CS); fn = fn(strncmp('a',fn,1));

end
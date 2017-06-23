function fn = cs_fieldnames(CS)
% CS FIELDNAMES - get antenna field names in Cross Spectra File
% fn = cs_fieldnames(CS)
% 
%     fn = {'antenna1Self'
%         'antenna2Self'
%         'antenna3Self'
%         'antenna12CrossSp'
%         'antenna13CrossSp'
%         'antenna23CrossSp'
%                           };


fn = fieldnames(CS); fn = fn(strncmp('a',fn,1));

end
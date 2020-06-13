function [fn,I,J] = cs_make_field_names(M)
% CS MAKE FIELD NAMES - create names for fields for non-seasonde CS
% [fn,I,J] = cs_make_field_names(M)
%
% Names for auto and cross spectra of the upper triangle 
%
% see also, cs_struct.m, cs_from_fft.m, cs_fieldnames.m, doa_column_index.m

% Copyright (C) 2017 Brian Emery



% just need the upper triangle 
[I,J] = find(triu(ones(M)));


% Zero pad if M > 9 (num2str(I,'%02.0f')_
if M > 9
    fmt = '%02.0f';
else
    fmt = '%1.0f';
end

fn = cell(length(I),1);

for n = 1:length(I)
        
        % make the new field name
        fn{n} = ['a' num2str(I(n),fmt) num2str(J(n),fmt) ];
        
end



end

function [desc,units,fmt] = ctf_descriptions(HdrNames,HdrValues)
% CTF_DESCRIPTIONS.M - get or save descriptions for CTF file keys
% [desc,units,fmt] = ctf_descriptions(HdrNames,HdrValues)
%
% This function is typically called by ctfReader.m (to populate the database
% as files are loaded), and ctfWriter.m to use info in the database.
%
% The 'database' is currently ctf_formats.mat, and saved on the function's
% matlab path in /codar/private/. The 'descriptions' are the strings associated
% with the 4 char keys, and their units.
%
% If necessary, this database can be recreated with code in ctf_formats.m,
% and or code in the test for ctfReader.m

% Copyright(C) 2012 Brian M. Emery
% Version 22 Feb 2012: 
% 	from ctf_formats.m

if nargin > 0
        
   % get descriptions from inputs
   [desc,units] = read_descriptions(HdrNames,HdrValues); 
      
   % add them to database
   save_ctf_formats(desc, 'desc' )
   save_ctf_formats(units,'units')
   
else
    
    % output descriptions and units ... and formats
    load('ctf_formats.mat','desc','units','fmt');    
    
end


end

function [D,U] = read_descriptions(HdrNames,HdrValues)
% READ DESCRIPTIONS
%
% parse header names and values to assign descriptions to each of the 4
% character keys
% 
% second output struct commonly contains unit descriptions

% get the field names to use from 'TableColumnTypes'
[~,~,val]=getNameValuePair('TableColumnTypes',HdrNames,HdrValues);
fn = cellstr(strparser(val{1}));

% Get the column names (possibly 2 rows here); after TableStart
ii = find(strcmp('TableStart',HdrNames));


% get, parse and pack colum data
D = pack_descriptions(fn,HdrNames{ii+1});

%keyboard

% GET UNITS
% sometimes this is where the units are located
% try for row 2 if it exists 
if numel(HdrNames) > ii+1
    
    U = pack_descriptions(fn, HdrNames{ii+2});
     
else
    
    % otherwise create the empty struct
    % U = cell2struct(cell(size(fn)),fn,1);
    for i = 1:numel(fn) 
        U.(fn{i}) = {''};
    end
   
end


end

function D = pack_descriptions(fn, HdrNames)

% strip off comment characters
HdrNames = regexprep(HdrNames,'%*','');

% get and parse colNames
colNames = cellstr(strparser(HdrNames)); 

% put the data in a struct
for i = 1:numel(fn)
    try
        D.(fn{i}) = colNames(i);
        
    catch E
        D.(fn{i}) = {''};
        
    end
end



end


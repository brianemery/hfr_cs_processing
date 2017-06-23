function fmt = ctf_formats( D, fn)
% CTF_FORMATS.M - get or save format strings for CTF file keys
% fmt = ctf_formats(fmt)
%
% Given a char input with somewhat regularly formatted columns, and and array
% of corresponding to 4 letter CTF column name keys, this function determines
% the corresponding string format values and saves them in a database. 
%
% Given no input and an output, this function will output the contents of
% the database in a struct with field names corresponding to 4 letter keys
% and values corresponding to the string formats. 
%
% This function is typically called by ctfReader.m (to populate the database
% as files are loaded), and ctfWriter.m to use info in the database.
%
% The 'database' is currently ctf_formats.mat, and saved on the matlab path
% 
% See code in test_case subfunction of ctfReader.m for adding NEW 4
% character key meta data by loading a loop file in the test data (a file
% named LOOP_new1_000000_000000.loop). 

% Copyright(C) 2012 Brian M. Emery
% Version 13 Feb 2012: 
% 	Combined this with read_formats subfunction so that it's all contained
%   Finalized outputs with nargin ==0 inputs

% TO DO
% Finish the 'output formats' case
% Write data: maybe just fprintf, as a .m? Currently going to dropbox .mat
%   file. (would be better for deployed code and hand editing)
% Case of same  keys with diff formats? build tests with lots of ctf files
%   - might allow for file type be included if diff #'s for same keys
% Put together html examples for sourceforge repositiory

% for i = 1:numel(fn), if isempty(units.(fn{i})), units.(fn{i}) ={''}; end

if nargin > 0
    
   % get formats from inputs
   fmt = read_formats( D, fn);
   
   % add them to database
   % save_formats(fmt)
   save_ctf_formats(fmt,'fmt')
   
else
    
    % output formats
    S = load('ctf_formats.mat');
    fmt = S.fmt; 
    
end


end

function fmt = read_formats( D, fn)
% READ FORMATS 
% totally experimental at this point - machine learning?
%
% doc sprintf for details re formatting
%
% The FORMAT string is of the form:  %<WIDTH>.<PREC><SPECIFIER>
%         <SPECIFIER> is required; <WIDTH> and <PREC> are optional.
%         <WIDTH> is the number of characters or digits to read.
%         <PREC> applies only to the family of %f specifiers, and specifies 
%         the number of digits to read to the right of the decimal point. 
 
% Why does D have sci notation like numbers in it?
% How would this fit into matlab's text reading code?
% This needs to get called by ctfReader, and then update a data base
% (subfunction?) called by ctfWriter.m
% test on loop, RDLm, STAT, etc 



% GET WIDTHS

% get location of first white space in char array
% s is array of cells
% - this will do the whole matrix
% s = regexp(cellstr(D),'\s+'); 
%
% s should be regular, but maybe we need to check?
% for now, just grab one row
s = regexp(cellstr(D(1,:)),'\s+');

% convert s into a matrix (or vector)
s = cell2mat(s);

% here're the field widths
% NEED to pad the leading blanks though (+3)
widths = diff([0 s+2 size(D,2)+1+2]); % one for each field name



% GET PRECISIONS
% parse numbers in each line into cells, then get the position of the '.'
% relative to the start of the cell, the compute the index of the '.' from
% the end - that's the precision

%  use textscan to parse each line into cells 
% (do for whole thing by looping I think ...)
% CTF files are aligned on the '.', so probably only need to do this for
% the first line ...
A = textscan(D(1,:),'%s');

% this gets the index of the . relative to the start of the number
k = strfind(A{:},'.');

% get sizes of each cell
cols = cellfun('size',A{:},2);

% find cells with no decimal in them, will eventually set the precisions = 0
% for these
rows = cellfun('size',k,1); 
k(rows == 0) = {NaN}; 

% make k a vector
k = [k{:}];

% compute precision field widths! (NaN if prec == 0)
prec = cols' - k;

% might be redundant, but make the NaN's zero
prec(isnan(prec)) = 0;


% create format strings (eg:'%8.1f')
for i = 1:numel(fn)
    fmt.(fn{i}) = ['%' num2str(widths(i)) '.' num2str(prec(i)) 'f'];
end



return

% NOTES


% find the locations of the .'s
k = strfind(D(1,:),'.');

% find the location of the first blank, by matching one or more occurances
% of white space
s = regexp(D(1,:),'\s+');


% save /Users/codar/Dropbox/ctfReader.mat


end

function custom_edit_code

% load /m_files/tools/codar/private/ctf_formats.mat

fmt.A1SN = '%9.2f';
fmt.A2SN = '%7.2f';
fmt.A3SN = '%7.2f';

fmt.SBG1 = '%10.2f';
fmt.DPRV = '%14.6f';
fmt.TIME = '%8.0f';
fmt.TRGB = '%8.3f';
fmt.TGRV = '%14.6f';

units.TRGB = {'cwN'};

fmt.A13M = '%10.5f';
fmt.A23M = '%10.5f';
fmt.AR3D = '%10.5f';

fmt.A13P = '%11.5f';
fmt.A23P = '%11.5f';
fmt.AR3P = '%8.1f';



desc.A13M = {'A13M'};
desc.A13P = {'A13P'};
desc.A23M = {'A23M'};
desc.A23P = {'A23P'};
desc.AR3D = {'A33M'};
desc.AR3P = {'A33P'};

units.A13M = {'Volts'};
units.A13P = {'Deg'};
units.A23M = {'Volts'};
units.A23P = {'Deg'};
units.AR3D = {'dBm'};
units.AR3P = {'Deg'};

desc.A2SN = {'A2SN'};
desc.A3NF = {'A3NF'};
desc.A3SN = {'A3SN'};

save /m_files/tools/codar/private/ctf_formats.mat

% USE THIS LOOP FILE TO POPULATE CTF_FORMATS
T = ctfReader(['/m_files/test_data/ctfReader/LOOP_new1_000000_000000.loop']);
T = ctfReader(['/m_files/test_data/ctfReader/LOOP_sni1_20120223_000000.loop']);



end
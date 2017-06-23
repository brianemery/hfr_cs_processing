function created = creation_info
% CREATION INFO - make structure with mfilename, path, and date meta data
% created = creation_info
% 
% Simplifies the documentation of structures and *.mat files

% Copyright (C) 2011 Brian Emery
%   March 29 2011

[st,ix]=dbstack('-completenames');  
[pth,nm,ext] = fileparts(st(end).file);
           

created.by = 'Brian Emery'; 
created.date = datestr(now,0);
created.with = [nm ext]; % calling mfile and directory
created.in = pth;


end

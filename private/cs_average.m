function AV = cs_average(CS)
% CS AVERAGE - compute mean of fields in multi element CS structs
% CSA = cs_average(CS)
%
% Uses all the data in the multi-element input ... robust to empties
% 
% For converting CSQ to CSS, or CSS to CSA for example
%
% Given input in dB, with phases resulting from a complex number converted
% to dB, this will compute the average of both, such that the result can be
% converted back to volts^2. "This is one of the simplest examples of 
% statistics of non-Euclidean spaces." From: 
% https://en.wikipedia.org/wiki/Mean_of_circular_quantities
%
% SEE ALSO
% cs_dbm2volts.m and related

% Copyright (C) 2017 Brian Emery
% 
% Originally From the COS simulation code, but now much different

% TO DO
% rm filename 


% check for test case
if strcmp('--t',CS), test_case, return, end


% get field names
fn = cs_fieldnames(CS);

% create output
AV = CS(1);


for i = 1:numel(fn)
    
    % need this to get the size
    X = cat(3,CS.(fn{i}));
    
    % compute average
    AV.(fn{i}) = sum(X,3)./size(X,3);
end


% run on phases to if we have them
if isfield(CS,'Phases')
    
    % make a new Phases struct that is multi-elemnent
    Ph = [CS.Phases];   
    
    for i = 1:numel(fn)
        
        % need this to get the size
        X = cat(3,Ph.(fn{i}));
        
        % compute average of circular quantities
        AV.Phases.(fn{i}) = atan2(sum(sin(X),3),sum(cos(X),3));
        
    end
end


end


function test_case
% TEST CASE
% 
% test case directory: /m_files/test_data/
%
% Get data together for this
% wd = '/projects/soo_apm/data/Data_COP1/for_figure/';
% 
% flist = get_file_list(wd,'CSQ*');
% 
% for i = 1:numel(flist)
%     CS(i) = ReadCS(flist{i});
% end
% 
% save /m_files/test_data/cs_ship_rm.mat CS

% TEST for use of Phase info for dB input

% Get data
load /m_files/test_data/cs_ship_rm.mat CS

CS = cs_volts2dbm(CS);

tic
CSA = cs_average(CS);
toc

keyboard


% FUNCTIONAL TEST - Compare with previous

% % NOTES
% tic, CSA = cs_average(CS); toc
% Elapsed time is 0.013143 seconds.


% Get data
load /m_files/test_data/cs_ship_rm.mat CS


[CS(1:5).ProcessingSteps] = deal({''});
[CS(1:5).Units] = deal({''});

% empty out the center and see if averaging can handle it
CS(3) = cs_struct(1);



% COS VERSION
% Elapsed time is 0.053224 seconds.
% tic
% CS = cs_average(CS,numel(CS));
% toc
% Elapsed time is 0.018301 seconds.
% tic
% CS = cs_average(CS,numel(CS));
% toc
% Elapsed time is 0.017488 seconds.


% NEW VERSION
% Elapsed time is 0.030244 seconds.
% cs_average('--t')
% Elapsed time is 0.024711 seconds.
% 
% ... similar performance: check
tic
CSA = cs_average(CS);
toc

H = cs_plot(CS(1),12); hold on,

LS = line_style_groups;

h = cs_plot(CSA,12,LS(3),H);

h = cs_plot(CS(5),12,LS(2),h);


keyboard

% EXTRA
% test the 3d mean of the phases
x = [-175 175 120 145 -177]; % arithmetic mean is ~77 which is not right!
for i = 1:5
    X(1:10,1:10,i) = x(i)*pi/180;
end

X_ = atan2(sum(sin(X),3),sum(cos(X),3)); % gets about 162 deg


end

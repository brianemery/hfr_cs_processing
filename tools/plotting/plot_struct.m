function [hx,h] = plot_struct(S,xstr,fn,LS,hx)
% PLOT STRUCT - plotting template for structures
% [hx,hdl] = plot_struct(S,xstr,fn,LS)
%
% Makes line plots of data in structure, such as timeseries, etc.
%
% INPUTS
% S    - data structure
% xstr - x axis field name string
% fn   - cell of field names to plot on y axis (multi element for subplots)
% LS   - line style structure(s), see line_style_groups.m
% hx   - optional axes handles to previous use of plot_struct
%
% If units and a LABS fields are found, these are used to label the axes
% ** need to specify these and standardize
%
% EXAMPLE
% md = '/projects/drifter_experiment_ssd/data/seasonde_day6/'
% apm_file = [md 'RadialConfigs/MeasPattern.txt'];
% APM = load_pattern_file(apm_file);
% M = interp_apm(APM,min(APM.BEAR):0.1:max(APM.BEAR));
% 
% LS = get_linestyles(8);
% [hx,hdl] = plot_struct(APM,'BEAR',{'A13I','A13R','A23I','A23R'},LS(1));
% hold on
% LS(2).Marker = '.'
% LS(2).LineStyle = 'none';
% [hx,hdl] = plot_struct(M,'BEAR',{'A13I','A13R','A23I','A23R'},LS(2),hx);
% 
% figure
% [hx,hdl] = plot_struct(APM,'BEAR',{'A13M','A13P','A23M','A23P'},LS(1));
% hold on
% [hx,hdl] = plot_struct(M,'BEAR',{'A13M','A13P','A23M','A23P'},LS(2),hx);
% 
%
% see also: saveFig.m

% Copyright (C) 2011 Brian M. Emery
% 10 Nov 2011

% TO DO
% add legends?
% modify figure name
% test cases
% generate the code to make the plot 'by hand'

% check inputs, 
if ischar(fn)
    fn = cellstr(fn);
end


% get line styles if not given
if nargin <4, 
    LS = get_linestyles(numel(fn));
end

if numel(fn) > numel(LS)
    [LS(1:numel(fn))] = deal(LS);
end


% MAIN PLOTTING 

if nargin < 5

    figure
    hx = makesubplots(numel(fn),1,.05,.05);

    % add letter
    subplot_add_letters(hx) %,0.035,0.85)

end

for i = 1:numel(fn)
        
    % plot data
    h{i} = plot(hx(i), S.(xstr), S.(fn{i}), LS(i)); hold on
    
    % label using field name (try units also)
    set_ylabel(hx(i),S,fn{i})
    
    % axes labels on bottom only
    if i < numel(fn)
       set(hx(i),'XTickLabel','')
    end
       
end


% set y limits automatically
set_ylims_auto(hx)

% set plot size and fonts, etc
publicationStandards




end

function set_ylabel(hx,S,fn)
% SET YLABEL
% try to use Units field, or LABS

% set default first
lab = fn;

% try several ways of storing units info (this is getting messy)
if isfield(S,'Units') && ~isfield(S,'LABS')
    try  lab = [fn ' ' S.Units.(fn)]; catch E,  end

elseif isfield(S,'LABS')
    
    % label storage standard
    try  lab = [fn ' ' S.LABS.(fn)]; catch E, end
    
    % label storage used in Sally project code
    try  lab = [fn ' ' S.LABS.(fn).u]; catch E, end   
    
    % label storage used by ctfReader.m and related
    try  lab = [S.LABS.(fn){:} '' S.Units.(fn){:} ]; catch E, end    
    
    
end

ylabel(hx,lab)

end

function LS = get_linestyles(n)
% GET LINESTYLES
%
% This might be the best way to do this, use the field name as the string
% THIS IS HOW WE DO IT
% takes advantage of structure inputs to plot!!!

% LINE STYLE SETTINGS

% % defaults
% LS(1).Color = 'k';
% LS(1).Marker = '.';
% % LS(1).MarkerFaceColor = 'w';
% % LS(1).MarkerEdgeColor = 'k';  
% LS(1).LineStyle = '-';
% % LS(1).MarkerSize = 6;
% 
% LS(1).Color = 'k';
% LS(1).Marker = 'o';
% LS(1).MarkerFaceColor = [.8 .8 .8]; %'w';
% LS(1).MarkerEdgeColor = 'k';
% LS(1).LineStyle = 'none';
% LS(1).MarkerSize = 4;
% 
% [LS(1:n)] = deal(LS(1));
% 
% LS(2).Color = 'r';
% LS(2).Marker = '.';
% LS(2).MarkerFaceColor = 'r';
% LS(2).MarkerEdgeColor = 'r';

c = line_colors;
ls = {'*','s','o','^','d','+','x'};


for i = 1:min([n numel(ls)])
    
    LS(i).Color = c(i,:);
    LS(i).Marker = ls{i};
    LS(i).MarkerSize = 6;
    LS(i).MarkerFaceColor = c(i,:);
    LS(i).MarkerEdgeColor = c(i,:);
    LS(i).LineStyle = '-';
    LS(1).LineWidth = 0.5000;

end



end

function LABS = get_y_info
% EXAMPLE 
 
% FIGURE 2

LABS.LE_0900.y = 'LE ';%_0_9_0_0';
LABS.LE_0900.u = '(W m^-^2)';
LABS.LE_0900.y_lim = [-350 0];

LABS.SE_0900.y = 'SE ';%_0_9_0_0';
LABS.SE_0900.u = '(W m^-^2)';
LABS.SE_0900.y_lim = [-40 10];

LABS.LWNet_0900.y = 'LWNet ';%_0_9_0_0';
LABS.LWNet_0900.u = '(W m^-^2)';
LABS.LWNet_0900.y_lim = [-100 -10];


end


function h = pws_map(ax,swch)
% PWS_MAP - map of Prince William Sound
% h=pws_map(ax,swch)
%
% Max axis for present coast line is : 
% ax=[-148  -146    60    61];
%
% INPUTS:
% optional, axes, and 'm' to use m_map toolbox
%
% 
% Needs bathymetry ...
% m_elev('contour',[-3500:1000:-500],'edgecolor','b');

% ----------------------------------------------------------------------
% DEFINITIONS
% Data directory
coast_dir = '/projects/drifter_simulation/pws/data/' ;

% Site locations
% SiteSource:  1 KNOW   60.6878667 -146.5176167   75.00     1.490 Meas     255.0
% SiteSource:  2 SHEL   60.4442667 -146.6623333   75.00     1.490 Meas     302.0
%
% And  Other Site Info
INFO.SHEL.str = 'Shelter Bay';
INFO.SHEL.x   = 0.03;
INFO.SHEL.y   = 0;
INFO.SHEL.loc = [-146.6623333 60.4442667];

INFO.KNOW.str = 'Knowles Bay';
INFO.KNOW.x   = -0.075;
INFO.KNOW.y   = +0.03;
INFO.KNOW.loc = [-146.5176167 60.6878667];


% -----------------------------------------------------------------------
% SORT INPUTS
% -----------------------------------------------------------------------
% set defaults
%;%[-147-40/60 -146-20/60 60+15/60 60+50/60];
default_ax = [-148  -146    60    61]; %[-147.8580 -146.3082   60.1740   60.9489]; 
swch=0;
h=[];

if nargin==0
    ax=default_ax;
    swch=0;
elseif strcmp(ax,'m')
    ax=default_ax;
    swch=1; % use m_maps tools if true
elseif nargin==2
    swch=1;
end


% ----------------------------------------------------------------------
% MAKE MAP
% ----------------------------------------------------------------------


if swch
    % M_MAP VERSION
    make_pws_mmap(coast_dir, INFO, ax);
else
    % NON-M_MAP VERSION
    h=make_pws_map(coast_dir, INFO, ax);
end

end


% ----------------------------------------------------------------------
function make_pws_mmap(coast_dir, loc, ax)

%load in SS bathymetry data
%result of running ss_topo_proc;
%load /Users/carter/Desktop/philex_iop09/data/core_data/phil_1min;


orient portrait;
m_proj('mercator','lon',[ax(1) ax(2)],'lat',[ax(3) ax(4)]);

% % make a coastline file and save it for later use
% % must be updated if map boundaries change
% m_gshhs_h('save','coast_pws');

%fill in land
m_usercoast([coast_dir 'coast_pws.mat'],'patch',[0.75 0.75 0.75],'linewidth',0.7); %
hold on;

m_grid('contour','on','tickdir','out','box','fancy','fontsize',14,'linewidth',1,'linestyle','none'); %'--'

% add hf sites
add_sites(loc,1)

%m_grid('contour','on','tickdir','out','box','on', ...
       %'xtick',[122.5:0.5:126.5], ...
       %'xticklabels',['122.5';'123.0';'123.5';'124.0';'124.5';'125.0';'125.5';'126.0';'126.5'], ...
       %'ytick',[8 9 10 11], ...
       %'yticklabels',[' 8^oN';' 9^oN';'10^oN';'11^oN'], ...
%       'fontsize',14,'linewidth',2,'linestyle','--');
% xlabel('Longitude','Fontsize',14,'Fontweight','bold');
% ylabel('Latitude','Fontsize',14,'Fontweight','bold');

end

% ----------------------------------------------------------------------
function h = make_pws_map(coast_dir, loc, ax)
% NON-M_MAP version
% 
% Good for quick and dirty plots but may have distortion in
% x vs y distances. Also allows alot of axes adjustments without too much
% trouble

load([coast_dir 'coast_pws.mat']);

hold on
h=plot(ncst(:,1),ncst(:,2),'-k');
axis(ax);

% % Add sites
% sites={'KNOW','SHEL'};
% for i = 1:numel(sites)
%    hs=plot(loc.(sites{i})(1),loc.(sites{i})(2),'k^');
%    text(loc.(sites{i})(1)+.05,loc.(sites{i})(2),sites{i})
%    set(hs,'MarkerFaceColor','y')
% end
add_sites(loc)

% Get mercator scaling factor the area and scale the map and area 
% (thanks to Mike Cook)
limits = axis;
[ax,sc]=mercat(limits(1:2),limits(3:4));
set(gca, 'DataAspectRatio', [1,sc,1],'PlotBoxAspectRatio',[1,1/ax,1]);

xlabel('Longitude','Fontsize',14,'Fontweight','bold');
ylabel('Latitude','Fontsize',14,'Fontweight','bold');

end

% ----------------------------------------------------------------------
function add_sites(INFO,mmap)
% ADD SITES


% Determine if using mmap tools
if nargin < 2
    mmap ='';
else
    mmap = 'm_';
end

% Add HF sites to plot
sites = fieldnames(INFO);

for i = 1:numel(sites)
    
    % use subfunction to clean up code
    plot_and_label_site(INFO.(sites{i}),mmap)
    
end


end

% ----------------------------------------------------------------------
function plot_and_label_site(S,mmap)
% PLOT AND LABEL SITE

% Marker Settings
LS.Color = 'k';
LS.Marker = '^';
LS.MarkerFaceColor = 'y';
LS.MarkerEdgeColor = 'k';  
LS.LineStyle = 'none';
LS.LineWidth = 0.9;
LS.MarkerSize = 10;


% Append 'm_' if necessary
eval([mmap 'plot(S.loc(1)     , S.loc(2)      , LS);'])
eval([mmap 'text(S.loc(1)+ S.x, S.loc(2) + S.y, S.str)'])

               
end
    


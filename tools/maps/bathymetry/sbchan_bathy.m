function h = sbchan_bathy(dpth)
% SBCHAN_BATHY.M
% Called by sbchan_map, now plots data for a large section of the so cal
% bight.

% Brian Emery
% originaly mid to late 1990's, 7dec09 updated

% Make a fading color map. This could be better
% cmap=[ones(length(d2)+1,1) (0:(1./length(d2)):1)' (0:(1./length(d2)):1)']; 

hold on 

% --------------------------------------------------------
%  DEFINITIONS
%----------------------------------------------------------
% Other possible settings for the color, etc of the bathy lines
%opt='''k'',''LineWidth'',0.5,''Color'',[.75 .75 .75]';
%opt='''k'',''LineWidth'',0.5';
%opt='''k'',''LineWidth'',0.5,''Color'',[1.0000    0.8235    0.8235]';

% DEPTHS TO PLOT
% define the depths to plot (Melanie used "20 m, 50 m, 100 m, and then every 100 m.") 
if nargin < 1
    dpth=[50 100:100:1000 1200:200:3600]';
end

% BATHY LINE COLOR
% For example, [0 0 0] is black, [1 1 1] is white,
%     [1 0 0] is pure red, [.5 .5 .5] is gray
%lineColor='k';
%lineColor = lineColor*ones(1,3);

% Make a fading color map like Melanie's
lineColor = zeros(length(dpth),3);
lineColor(dpth<=2500,:) = interp1([50 2500]',[0.75 0.75 0.75; 0 0 0],dpth(dpth<=2500,:)); %[0.8481 0.8481 0.8481; 0 0 0],dpth);

% all these should be in a folder on the path (eg )
fileToLoad={'sbcsmb_quick.mat'; ...
    'la_quick.mat'; ...
    'mugu_quick.mat'; ...
    'outer_quick.mat'; ...
    'soChanIs_quick.mat';};

% --------------------------------------------------------
%  LOAD AND PLOT 
%----------------------------------------------------------

for i=1:length(fileToLoad)
    h{i} = bathy_plotter(fileToLoad{i},dpth,lineColor);
end

% % Add Labels to Contours
% label_contours(lineColor);
% 

end

% --------------------------------------------------------
function h = bathy_plotter(fileToLoad,dpth,lineColor)
% BATHY PLOTTER
% subfunction to do all the work

load(fileToLoad)

% Keep handles
h=NaN(size(dpth));

% convert to strings
d=cellstr(num2str(dpth));

for ii=1:length(d)
   eval(['h(ii)=plot(d' strtrim(d{ii}) '(:,1),d' strtrim(d{ii}) '(:,2),''color'',[' num2str(lineColor(ii,:)) ']);'],'s=1;')
end

h = h(~isnan(h));

% % Add some text 
% ht=text(-119.8358,33.8495,'bathymetry in meters','fontsize',6,'color',lineColor(1,:));

% Set clipping for all bathy objects
set(h,'Clipping','on')

end
% --------------------------------------------------------
function label_contours(lineColor)
% LABEL CONTOURS
% this is all totaly hard wired, hence the seperate function
% All these postions obtained with ginput

% Brian Emery

% label bathy contours
% 'labs' is obtained using ginput
dlab =num2str([50 100:100:600 50 100:100:600 ...
    700 800 1000 1400 1800 2000 2200]');

% Lab positions for d<700 (twice)
labxy=[-120.6091   34.5172
 -120.5171   34.4197;
 -120.5133   34.3872;
 -120.5152   34.3624;
 -120.5208   34.3299;
 -120.5603   34.2757;
 -120.6279   34.2386;
 -119.6953   34.0444;
 -119.7461   34.0912;
 -119.7791   34.1159;
 -119.7791   34.1356;
 -119.8270   34.2120;
 -119.9258   34.2836;
 -120.0275   34.2367;];

% Lab postions d>600 and d<2300
labxy=[labxy;
 -120.6975   34.2030;
 -120.7426   34.0594;
 -120.7407   34.0285;
 -120.7407   33.9992;
 -120.7445   33.9729;
 -120.7632   33.8949;
 -120.8046   33.8602;];

% put the Labels on the contours
for i = 1:length(dlab)
    ht=text(labxy(i,1),labxy(i,2),strtrim(dlab(i,:)),'fontsize',6,'color',lineColor(i,:));
end

% % label bathy contours
% % 'labs' is obtained using ginput
% d=[d' d']';
%  %labs=[-120.7280   34.8359; %these put labs in sm basin
%  %-120.7999   34.8334;
%  %-120.8807   34.8285;
%  %-120.9585   34.8285;
%  %-121.0184   34.8014;
%  %-121.0872   34.7742;
%  %-121.0812   34.6731;
%  labs=[-120.6091   34.5172
%  -120.5171   34.4197;
%  -120.5133   34.3872;
%  -120.5152   34.3624;
%  -120.5208   34.3299;
%  -120.5603   34.2757;
%  -120.6279   34.2386;
%  -119.6953   34.0444;
%  -119.7461   34.0912;
%  -119.7791   34.1159;
%  -119.7791   34.1356;
%  -119.8270   34.2120;
%  -119.9258   34.2836;
%  -120.0275   34.2367;];
% 
% % put the Labels on the contours
% for i=1:length(labs)
%    ht=text(labs(i,1),labs(i,2),d(i),'fontsize',6,'color','r');
%       set(ht,'Clipping','on')
% end

% % add deep countours. The data is in the data file already, but this
% % can be commented out if required
% 
% d2={'700';'800';'900';'1000';'1200';'1400';'1600';'1800'; ...
%       '2000';'2200';'2400';'2600';'2800';'3000';'3200';'3400';'3600';};
% 
% % Make a fading color map. This could be better
% % cmap=[ones(length(d2)+1,1) (0:(1./length(d2)):1)' (0:(1./length(d2)):1)']; 
% opt='''k'',''LineWidth'',0.5,''Color'',[1.0000    0.8235    0.8235]';
% 
% for ii=1:length(d)
%    % Use the top command with fading colormap
%    %eval(['plot(d' d2{ii} '(:,1),d' d2{ii} '(:,2),''Color'',[' num2str(cmap(ii,:)) '])'])
%    eval(['plot(d' d2{ii} '(:,1),d' d2{ii} '(:,2),' opt ')'])
%    %eval(['plot(d' d2{ii} '(:,1),d' d2{ii} '(:,2),-' d2{ii} '*ones(size(d' d2{ii} '(:,1)))/1000,' opt ')'])
% end
% 
% % % add this 
% ht=text(-119.8358,33.8495,'bathymetry in meters','fontsize',6,'color','r');
% set(ht,'Clipping','on')
% 
% % label bathy contours
% % 'labs' is obtained using ginput
% labs=[-120.6975   34.2030;
%  -120.7426   34.0594;
%  -120.7407   34.0285;
%  -120.7407   33.9992;
%  -120.7445   33.9729;
%  -120.7632   33.8949;
%  -120.8046   33.8602;];
% 
% % put the Labels on the contours
% d3={'700';'800';'1000';'1400';'1800'; '2000'; '2200';};
% for i=1:length(labs)
%    ht=text(labs(i,1),labs(i,2),d3(i),'fontsize',6,'color',[1.0000    0.8235    0.8235]);
%    set(ht,'Clipping','on')
% end

% Set clipping for all bathy objects
set([ht(:);],'Clipping','on')


end




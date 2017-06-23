function percent=total_coverage_map(U,serialtime,AngDiff,ang_dif,gridd,DIR)
% TOTAL COVERAGE MAP.M
% percent = total_coverage_map(U,serialtime,AngDiff,ang_dif,gridd,DIR)
% Computes and plots percent coverage given total vector data.
%
% DIR is optional and needs some documentation (real time processing)
%
% see also using plotData from HFRP tools, eg:
% handles = plotData( TUV, 'perc','plot'); does a coverage map
%
% TUV usage:
% p = total_coverage_map(TUV.U,TUV.TimeStamp,[],[],TUV.LonLat);

% 19Feb99 Brian Emery
% 30Oct03 Improved to use less memory

if nargin<6, DIR=[];, end

if ~isempty(AngDiff) && ~isempty(ang_dif)
    % First sort out the points that dont satisfy ang_dif
    U(AngDiff<=ang_dif(1) | AngDiff>=ang_dif(2))=NaN; % these take up all the time
    clear AngDiff ang_dif
end

% set the NaN's in U to zeros and append it to the coverage matrix
coverage=~isnan(U); clear U
coverage=+coverage; % make it a double array

% compute percent from the coverage matrix
[num_gridpts,hours]=size(coverage);
percent=(sum(coverage,2)./hours)*100; % this takes a lot of time too.
clear coverage

% blank the zero values
percent(percent==0)=NaN; 

% plot
% sbc, hold on, codar_axis
cdot2d(gridd(:,1),gridd(:,2), percent,10,'jet',[0 100]);  % 15 for dots
set(gcf,'units','pixels','position',[199   179   871   630])

% if ~isempty(DIR)
%     hleg=sbchan_map_auto(DIR); 
% else
%     sbc
% end

% save the title info as a string
timestr=['Percent Coverage']; % for ' datestr(serialtime(1),3) ' ' datestr(serialtime(1),10) ];
   
title(timestr);
xlabel('Longitude');
ylabel('Latitude');

% more stuff from plotrad2tot.m
set(gca,'FontSize',12,'FontWeight','bold');
set(get(gca,'Xlabel'),'FontSize',12,'FontWeight','bold');
set(get(gca,'Ylabel'),'FontSize',12,'FontWeight','bold');
set(get(gca,'Title'),'FontSize',12,'FontWeight','bold');
set(findobj('type','text'),'clipping','on');
%text(-119.8, 34.6,['Angdif= ' num2str(ang_dif)]); 

return

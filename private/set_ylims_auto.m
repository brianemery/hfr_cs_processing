function set_ylims_auto(hx)
% SET YLIMS AUTO - set figure y limits based on data
% set_ylims_auto(hx)
%
% Given axes handle(s) sets the y limits of the plot based on the data
% contained on the figure.

% Copyright(C) 2011 Brian Emery
% 11 Nov 2010
% 11 Nov 2011 - refactored, generalized

% TO DO
% - build test case, with multiple lines on the plot
% - option for subplots to have same (max) ylims
% see vic_climate_fluxes_noaa_summary_figure for a later modification

for i = 1:numel(hx)
    
    try
        
%         % get all x and y data that went into the plot
%         S = get(get(hx(i),'Children'));
%         
%         % concatenate it
%         x = [S.XData];
%         y = [S.YData];
        
        % get all x and y data that went into the plot
        [x,y] = get_data(hx(i));

        
        % get current xlimits
        xlims = get(hx(i),'xlim');
        
        % find y values of data between x limits (visible data)
        ix = find( x > xlims(1) & x < xlims(2) );
        
        % limit y data
        y = y(ix);
        
        % add 10% to final plot y limits
        dy = 0.05*diff([min(y) max(y)]);
        
        % try to set the limits
        ylims = [min(y)-dy  max(y)+dy];
        
        if all(~isnan(ylims))
            set(hx(i),'ylim',ylims)
        end
        
    catch err
        
        
        % display cause of any errors
        %keyboard
        disp([mfilename,':', err.message])
        
    end
    
    
end

% OLD WAY - pretty cool
% 
% set_ylims(hx(1), {D.temp})
%
% function set_ylims(hx,y)
% % SET YLIMS - dynamically set ylimits
% %
% %
% % get ymax (either 21 or max of stuff in D, +1)
% % This is just a fancy way to concat D.temp outputs
% 
% % S = get(get(hx,'Children'));
% % S(1)
% 
% try
%     ymax = max(cellfun(@max,y))+0.75;
%     ymin = min(cellfun(@min,y))-0.75;
%     
% catch % R2007 format
%     disp('ylim setting error')
%     ymax = 20; ymin = 10;
% end
% 
% set(hx,'ylim',[ymin ymax])
% 
% end
% 
% 


end

function [x,y] = get_data(hx)
% get all x and y data that went into the plot

% STUCK
% hard to separate the legend objects from the data objects ... needs more
% work

% get handles of all figure objects
hdls = get(hx,'Children');

% get just the data handles
str = get(hdls,'Type');

% get the struct containing object data
S = get(hdls(strcmp('line',str)));

        
% concatenate it
x = [S.XData];
y = [S.YData];
        
end


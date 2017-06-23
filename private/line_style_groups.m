function LS = line_style_groups(n)
% LINE STYLE GROUPS - get favorite line style groups
% LS = line_style_groups(n)
%
% Given an array of plot object handles, and a group name (see list below),
% this sets the color, linestyle, marker, and marker face color using some
% of my favorite collections.
%
% ... or something lijke that, just an aggregation of the above ...
%
% EXAMPLE
% LS = line_style_groups;
% [LS(1:2).LineStyle] = deal('-');
%
% NOTE
% new matlab default color order: (see color
%     [    0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840]


% TO DO
% more flexibility, customizeations?
%



% Copyright (C) 2010 Brian M. Emery


% THIS IS HOW WE DO IT
% takes advantage of structure inputs to plot!!!

if nargin < 1
    n = 2;
end



% LINE STYLE SETTINGS

% raw apm data
LS(1).Color = [.8 .8 .8];
LS(1).Marker = '.';
LS(1).MarkerSize = 6;
LS(1).MarkerFaceColor = 'none';
LS(1).MarkerEdgeColor = [.8 .8 .8];
LS(1).LineStyle = 'none';
LS(1).LineWidth = 0.5000;


% binned apm data
LS(2).Color = 'k';
LS(2).Marker = 'o';
LS(2).MarkerSize = 6;
LS(2).MarkerFaceColor = 'k';
LS(2).MarkerEdgeColor = 'k';
LS(2).LineStyle = '-';
LS(2).LineWidth = 0.5000;


% measured apm data
LS(3).Color = 'k';
LS(3).Marker = '.';
LS(3).MarkerSize = 6;
LS(3).MarkerFaceColor = 'k';
LS(3).MarkerEdgeColor = 'k';
LS(3).LineStyle = 'none';
LS(3).LineWidth = 0.5000;






return
% LINE STYLE SETTINGS


LS(1).Color = 'k';
LS(1).Marker = 'o';
LS(1).MarkerFaceColor = 'w';
LS(1).MarkerEdgeColor = 'k';  
LS(1).LineStyle = 'none';
LS(1).MarkerSize = 8;

% defaults
LS(2:n) = LS(1);

% [LS(1:n).LineStyle] = deal('-');
% [LS.LineWidth] = deal(1.5);
% [LS.MarkerSize] = deal(8);


% Custom for LS(2)
LS(2).Color = 'g';
LS(2).Marker = 'd';
LS(2).MarkerFaceColor = 'g';
LS(2).MarkerEdgeColor = 'g';

return


LS(1).Color = 'r';
LS(1).Marker = '.';
LS(1).LineStyle = '-';

LS(2).Color = 'b';
LS(2).Marker = '.';
LS(2).LineStyle = '-';

LS(3).Color = 'g';
LS(3).Marker = '.';
LS(3).LineStyle = '-';





% DEFINE LS (LineStyles) structure

% raw apm data
LS.A.Color = [.8 .8 .8];
LS.A.LineStyle = '.';
LS.A.LineWidth = 0.5000;
% LS.A.Marker = 'none';
% LS.A.MarkerSize = 6;
% LS.A.MarkerFaceColor = 'none';

% binned apm data
LS.B.Color = 'k';
LS.B.LineStyle = '-';
LS.B.LineWidth = 0.5000;
LS.B.Marker = 'o';
LS.B.MarkerSize = 6;
LS.B.MarkerFaceColor = 'k';

% measured apm data
LS.M.Color = 'k';
LS.M.LineStyle = ':';
LS.M.LineWidth = 0.5000;
LS.M.Marker = 'none';
% LS.M.MarkerSize = 6;
% LS.M.MarkerFaceColor = 'k';

% Store it
LABS.LS = LS;

% Set Typical Yaxis settings
vars = {'A13M','A23M','A13R','A23R','A13I','A23I'};
for i = 1:numel(vars)
    LABS.(vars{i})(1).ylims = [-1 1];
end
vars = {'A13P','A23P'};
for i = 1:numel(vars)
    LABS.(vars{i})(1).ylims = [-185 185];
end





return


% LINE STYLES
% each row a sites color, lineStyle, Marker, face color, marker edge color
% ... there must be a better way, hmmmm... maybe some way of taking
% advantage of the struct inputs ...
clmf = {'k'          ,'-'   ,'o'    ,'k'           ,'k'           ; ...
        'k'          ,':'   ,'o'    ,'w'           ,'k'           ; ...
        [.5 .5 .5]   ,':'   ,'v'    ,'w'           ,[.5 .5 .5]    ; ...
        [.75 .75 .75],':'   ,'s'    ,'w'           ,[.75 .75 .75] ; ...
        [.75 .75 .75],'-'   ,'s'    ,[.75 .75 .75] ,'k'           ; ...
        [.5 .5 .5]   ,'-'   ,'v'    ,[.5 .5 .5]    ,'k'           ; } ;



for i = 1:6
    h(i) = plot(MET.(xf)(idx,i),MET.(yf)(idx,i)); hold on
    %set(h(i),'color',clr{i},'linestyle',lns{i},'marker',mkr{i},'markersize',ms,'markerfacecolor',mfc{i});
    set(h(i),       'color',clmf{i,1}, ...
                'linestyle',clmf{i,2}, ...
                   'marker',clmf{i,3}, ...
          'markerfacecolor',clmf{i,4}, ...
          'markeredgecolor',clmf{i,5}, ...
               'markersize',ms);
end





set(h(i),'color',clr{i},'linestyle',lns{i},'marker',mkr{i},'markersize',ms,'markerfacecolor',mfc{i});

% Set line color, style, marker type and color (for 6 sites)
% order: {'entebbe','jinja','kisumu','musoma','mwanza','bukoba'};
clr = {'b'; 'b'; 'r'; 'g'; 'r'; 'g';};
lns = {'-'; ':'; ':'; ':'; '-'; '-';};
mkr = {'o'; 'o'; 'v'; 's'; 'v'; 's';};
mfc = {'b'; 'w'; 'w'; 'w'; 'r'; 'g';};






% WHAT CO LIKES BEST
% real
set(h.real,'LineStyle','-' ,'Color','k','Marker','o','MarkerFaceColor','w')

% 4 models: sim, noturb, fixed, propto
set(h.sim,'LineStyle','-' ,'Color','k') 
set(h.noturb,'LineStyle','-' ,'Color','k','Marker','s')         
set(h.fxdstd,'LineStyle','--' ,'Color','k') 
set(h.std_100m,'LineStyle','-.' ,'Color','k')                                      

% com
set(h.com,'LineStyle','-' ,'Color',[.6 .6 .6])     

% keeps the legend just small enough
set(legh,'FontSize',13)

% WHAT I LIKE BEST
% % real
% set(h(1),'LineStyle','-' ,'Color','k')                                           
% 
% % 4 models: sim, noturb, fixed, propto
% set(h(2),'LineStyle','-' ,'Color','k','Marker','o','MarkerFaceColor','w')         
% set(h(4),'LineStyle','-' ,'Color','k','Marker','s')         
% set(h(5),'LineStyle','--' ,'Color',[.8 .8 .8])
% set(h(6),'LineStyle','-' ,'Color',[.6 .6 .6])                                     
% 
% % com
% set(h(3),'LineStyle','-' ,'Color',[.8 .8 .8],'Marker','s')                                     




end

% alternate method

function LS = get_linestyles(str)
% GET LINESTYLES
%
% This might be the best way to do this, use the field name as the string
% input

switch str
    
    case 'R5'
        LS.Color = 'k'; %[.6 .6 .6];
        LS.MarkerFaceColor = [.6 .6 .6];
        LS.Marker = 'o';
        LS.MarkerSize = 6;
        LS.LineWidth = 0.7;
        LS.LineStyle = 'none';
        
    case 'RS'
        LS.Color = 'r';
        LS.MarkerFaceColor = 'r';
        LS.LineStyle = 'none';
        LS.Marker = 'o';
        LS.MarkerSize = 3;     
        
    case 'R1'
        LS.Color = 'c';
        LS.MarkerFaceColor = [.6 .6 .6];
        LS.LineStyle = 'none';
        LS.Marker = 'o';
        LS.MarkerSize = 3;     
        
    case 'DRFT'
        LS.Color = 'k';
        LS.LineStyle = 'none';
        LS.Marker = 'o';
        LS.MarkerFaceColor = 'k';
        LS.MarkerSize = 2;
        LS.LineWidth = 1;

    otherwise
        keyboard
end

end


function LS = line_style_groups_eg
% LINE STYLE GROUPS - get favorite line style groups

% LINE STYLE SETTINGS


LS(1).Color = 'r';
LS(1).Marker = 'o';
LS(1).MarkerFaceColor = 'r';
LS(1).MarkerEdgeColor = 'r';  
LS(1).LineStyle = '-';
LS(1).MarkerSize = 8;

% defaults
LS(2:4) = LS(1);


LS(2).Color = 'b';
LS(2).MarkerFaceColor = 'b';
LS(2).MarkerEdgeColor = 'b';  

LS(3).Color = 'm';
LS(3).MarkerFaceColor = 'm';
LS(3).MarkerEdgeColor = 'm';  

LS(4).Color = 'g';
LS(4).MarkerFaceColor = 'g';
LS(4).MarkerEdgeColor = 'g';  


end
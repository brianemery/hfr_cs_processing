function c = line_colors
% LINE COLORS - get default line colors
% c = line_colors .... or just use line.m from TMW
%
% c  =      [    0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840]
%
% blue, orange, yellow, purple, green, light blue, maroon
%
% SEE ALSO
% parula, eg:
% clr =  colormap_parula(10);
% figure, for i=1:10, plot([0 1],[i i],'-','Color',clr(i,:),'LineWidth',8), hold on, end


c =     [    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


end
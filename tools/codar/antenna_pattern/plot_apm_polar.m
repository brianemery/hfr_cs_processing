function H = plot_apm_polar(APM,LS,HH)
% PLOT APM POLAR - polar antenna patter plot
% h = plot_apm_polar(APM,LSH)
% Make Plots like Cross Loop Patterner
%
% INPUT
% Standard APM structure (minimum required fields are BEAR A13M A23M, see
% apm_struct.m)
%
% Optionally include Handle Struct (from previous run of this function) to
% enable adding second APM to plot
%
% EXAMPLE
% % Load and plot
% I = load_pattern_file('/m_files/test_data/make_ideal_pattern/IdealPattern.txt');
% H = plot_apm_polar(I); 
% 
% % Load, and over plot
% M = load_pattern_file('/m_files/test_data/make_ideal_pattern/MeasPattern.txt');
% LS = line_style_groups;
% H = plot_apm_polar(M,LS,H);


% Copyright (C) 2009-2010 Brian M. Emery
% June 2009
% Feb 2015 updates to make it more useful


if ~isfield(APM,'BEAR'), try APM.BEAR = APM.TRGB; catch, end, end

if nargin < 2
    LS = line_style_groups;
end


% FIRST PLOT POLAR

if nargin == 3
    figure(HH(1).fig)
else
    H(1).fig = figure;
end

H(1).h{1} = polar(APM.BEAR*pi/180,APM.A13M,LS(1)); hold on
H(1).h{2} = polar(APM.BEAR*pi/180,APM.A23M,LS(2));

try
    H(1).loop1 = polar([APM.loop1Brg APM.loop1Brg]*pi/180,[0 max([APM.A23M(:)' APM.A13M(:)'])],'y');
catch
end
view(90,-90)

%H(1).leg = legend_append([H(1).h{:}],{'A13M','A23M'});
H(1).leg = legend([H(1).h{:}],{'A13M','A23M'});

xlabel('Degrees CWN')



% SECOND PLOT OF PHASES

if nargin == 3
    figure(HH(2).fig)
else
    H(2).fig = figure;
end


H(2).h{1} = plot(APM.BEAR,APM.A13P,LS(1)); hold on
H(2).h{2} = plot(APM.BEAR,APM.A23P,LS(2));
xlabel('Bearing CWN')
ylabel('Phase (deg)')

H(2).leg = legend([H(2).h{:}],{'A13M','A23M'});




end


function LS = line_style_groups
% LINE STYLE GROUPS - get favorite line style groups

% LINE STYLE SETTINGS


LS(1).Color = 'r';
LS(1).Marker = 'o';
LS(1).MarkerFaceColor = 'r';
LS(1).MarkerEdgeColor = 'r';  
LS(1).LineStyle = '-';
LS(1).MarkerSize = 8;

% defaults
LS(2) = LS(1);


LS(2).Color = 'b';
LS(2).MarkerFaceColor = 'b';
LS(2).MarkerEdgeColor = 'b';  
% 
% LS(3).Color = 'm';
% LS(3).MarkerFaceColor = 'm';
% LS(3).MarkerEdgeColor = 'm';  
% 
% LS(4).Color = 'g';
% LS(4).MarkerFaceColor = 'g';
% LS(4).MarkerEdgeColor = 'g';  
% 

end
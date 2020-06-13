function H = plot_apm(A,LS,HH)
% PLOT APM - antenna patter plots like CrossLoopPatterner
% h = plot_apm(APM,LSH)
% Make Plots like Cross Loop Patterner - not just polar
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
% H = plot_apm(I); 
% 
% % Load, and over plot
% M = load_pattern_file('/m_files/test_data/make_ideal_pattern/MeasPattern.txt');
% LS = line_style_groups;
% H = plot_apm(M,LS,H);


% Copyright (C) 2009-2010 Brian M. Emery
% June 2009
% Feb 2015 updates to make it more useful


if ~isfield(A,'BEAR'), try A.BEAR = A.TRGB; catch, end, end

if nargin < 2
    LS = line_style_groups;
end


hx = makesubplots(2,2,.2,.2);

% phases
plot(hx(1), A.BEAR,A.A13P,LS(1)), hold on
plot(hx(1), A.BEAR,A.A23P,LS(2))

xlabel(hx(1),'Degrees (CWN)')
ylabel(hx(1),'Phases')

% amplitudes
plot(hx(2), A.BEAR,A.A13M,LS(1)), hold on
plot(hx(2), A.BEAR,A.A23M,LS(2))

% tangents
plot(hx(3), A.BEAR,tan(A.A13M/A.A23M),LS(1))

% real and imag 
plot(hx(4), A.BEAR,A.A13R,LS(1)), hold on
plot(hx(4), A.BEAR,A.A13I,LS(3))

plot(hx(4), A.BEAR,A.A23R,LS(2))
plot(hx(4), A.BEAR,A.A23I,LS(4))


% labels?

% handle managment

H = [];







return



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

H(1).leg = legend_append([H(1).h{:}],{'A13M','A23M'});

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

H(2).leg = legend_append([H(2).h{:}],{'A13M','A23M'});




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
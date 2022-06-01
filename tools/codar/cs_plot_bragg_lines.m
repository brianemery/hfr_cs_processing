function hx = cs_plot_bragg_lines(CS,h,a)
% CS PLOT BRAGG LINES - add bragg lines to CS plot
% 
% USAGE
% hx = cs_plot_bragg_lines(CS,h)
%
% OPTIONALLY
% add the ylims to specify height of the Bragg ticks
%
% see bragg_frequency_notes.m


% Brian Emery 

if nargin < 2, h = gcf; end


% constants
g = 9.81; %m/s^2
c = 299792458;% m / s

% ftx=13.49e6; % transmit freq in hz (1/s)
ftx = CS.Header.freqMHz *1e6; % transmit freq in hz (1/s)

ktx = 2*pi*ftx/c; % transmitter wave number

% Using the doppler relation for reflected light:
df = (1/(2*pi))*sqrt(2*g*ktx);

if nargin < 3
    a = axis(h);
    
else
    a = [NaN NaN a];
end

hx = plot(h,[df -df; df -df],[a(3:4)' a(3:4)'],'-k');





end

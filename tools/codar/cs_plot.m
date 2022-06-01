function H = cs_plot(CS,rcellIdx,LS,HH,ix)
% CS PLOT - plot one cross spectra's range cell
% H = cs_plot(CS,rcellIdx,LS,H,ix)
%
% LS is either a line style or struct of line styles
% eg:  LS = line_style_groups;
%
% fourth and fifth inputs allows this to be used to overplot a selection of doppler
% indecies
%
% Similar to:
% plot(CS.freqs,volts2dbm(CS.antenna3Self))
%
% ASSUMPTIONS
% CS matricies are nfft x nRcell in dimension
%
% USEFUL EXAMPLE
% fn ='/m_files/test_data/cs_plot/CSS_cop1_06_05_30_2100.cs4';
% CS = ReadCS(fn);
% rcIdx = 3;
% fbinIdx = 202;
% 
% % Make the plot
% figure
% H = cs_plot(CS,rcIdx); hold on,
% 
% LS = line_style_groups;
% 
% h = cs_plot(CS,rcIdx,LS(1),H,160);


% Copyright (C) 2011 Brian Emery 
%
% Revised 2014 into subplots, more labeling, etc

% TO DO
% - need ability to add to plot, pass the axes handle for htis ...
% - find code for adding axis with radial velocities on it (Pws?)


% Optionally run test case
if strcmp(CS,'--t'), test_case, return, end


% Sort out input options
if nargin < 3, 
    LS = default_linestyles;
end

if nargin < 5
    ix = ':';
end

if isstruct(LS) && numel(LS) < 3
    [LS(1:3)] = deal(LS(1));
elseif ischar(LS)
    [ls(1:3)] = deal(cellstr(LS));
    LS = ls;
end


[freqs,Vrad,dv] = getDopplerVelocities(CS.Header);



% Convert to dBm
CS = cs_volts2dbm(CS);

% Define fields to plot
fn = {'antenna13CrossSp','antenna23CrossSp','antenna3Self'};

% and the labels
fnlab = {'A13','A23','A33'};



% MAKE FIGURES
if nargin < 4
    hx = makesubplots(3,1,.05,.05);
    
else
    for i = 1:3
        hx(i) = HH(i).ax;
    end
end

for i = 1:3
    
    H(i) = plot_spectra(hx(i),freqs(ix),CS.(fn{i})(ix,rcellIdx),LS(i),fnlab{i});
    
    cs_plot_bragg_lines(CS,hx(i))

end



% Return early if just adding to plot
if nargin > 3
    H = HH; 
    return
end


% ANNOTATIONS

set(hx(1),'XTickLabel','')
set(hx(2),'XTickLabel','')

% Compute range
Range = ( 1:size(CS.antenna1Self,2) )*CS.Header.distToFirstRangeCell;

% Get the filename sans path (as a cell)
fname = fileparts_name(CS.FileName);

% Add title and info
ht = title(hx(1),[fname{1} ': Range Cell ' num2str(rcellIdx) ' (' num2str(Range(rcellIdx)) 'km)']);

set(ht,'Interpreter','none')

% x label
hxlab = xlabel(hx(3),'Frequency (Hz)'); %Doppler Velocity (cm/s)');



% This one is same for each subplt
% set_ylims(CS,fn,hx)
% This one is different for each subplot
set_ylims_auto(hx)

set(hx,'xlim',[min(freqs) max(freqs)])
% % set axes?
% if nargin < 3
%     for i = 1:length(hx)
%         axis(hx(i),[min(freqs) max(freqs) -165 -100])
%     end
% end


% SET OUTPUTS
[H(1:3).xlab]  = deal(hxlab);
[H(1:3).title] = deal(ht);


end

function H = plot_spectra(hx,freqs,CSM,LS,fnlab)
% PLOT SPECTRA

% Hack to deal with struct vs char line style inputs
if iscell(LS), LS = LS{:}; end

% simple plot
H.h = plot(hx,freqs,CSM,LS);


% Annotations     
H.ylab = ylabel(hx,[fnlab '(dBm)']);


% Outputs
H.ax = hx;

end

function LS = default_linestyles


LS(1).Color = 'r';
LS(1).Marker = '.';
LS(1).LineStyle = '-';

LS(2).Color = 'b';
LS(2).Marker = '.';
LS(2).LineStyle = '-';

LS(3).Color = 'g';
LS(3).Marker = '.';
LS(3).LineStyle = '-';

[LS(1:3).LineWidth] = deal(1.5);


end

function set_ylims(CS,fn,hx)


for i = 1:numel(fn)
    
    mx(i) = max(max(CS.(fn{i})));
    mn(i) = min(min(CS.(fn{i})));

end

mx = max(mx);
mn = min(mn);

for i = 1:3
    set(hx(i),'ylim',[mn mx])
end



end


function csq_plot_old(Vrad,CS,rcellIdx,velIdx,shipRadVel,shipRng,peakIdx,nBar,titleStr)
% ESTIMATE APM - PLOT SPECTRA
% A subfunction specific to compute_apm_from_csq.m (above), for making check plots
% Make a plot showing spectra w/ship info.

% PLOT SPECTRA
h1=plot(Vrad,10*log10(CS.antenna3Self(:,rcellIdx)),'-b.'); hold on
h5=plot(Vrad,10*log10(CS.antenna13CrossSp(:,rcellIdx))-20,'-r.'); hold on
h6=plot(Vrad,10*log10(CS.antenna23CrossSp(:,rcellIdx))-40,'-g.'); hold on

% PLOT SHIP LOCATION PATCH
% add ship location based on radial velocities
add_ship_to_csa(shipRadVel,shipRng)

% PLOT SHIP SIGNAL DETECTION
% Plot the signal that is suspected to be ship
h2=plot(Vrad(velIdx),10*log10(CS.antenna3Self(velIdx,rcellIdx)),'m*'); % ship!
h7=plot(Vrad(velIdx),10*log10(CS.antenna13CrossSp(velIdx,rcellIdx))-20,'m*');
h8=plot(Vrad(velIdx),10*log10(CS.antenna23CrossSp(velIdx,rcellIdx))-40,'m*');

% OTHER
h3=plot(Vrad(peakIdx{rcellIdx}),10*log10(CS.antenna3Self(peakIdx{rcellIdx},rcellIdx)),'y.');
ht=title(titleStr); set(ht,'interpreter','none')
a=axis;
hn=plot(a(1:2),[nBar nBar],'b');
xlabel('Doppler Radial Velocity (cm/s)')
legend([h5 h6 h1 h2 h3 hn], 'spectra 13 (-10db)','spectra 23 (-20db)', ...
    'spectra 33','ship?','1st order','noise level','Location','Best')

end

function hp = add_ship_to_csa(shipRadVel)
% ADD_SHIP_TO_CSA
% add_ship_to_csa(shipRadVel,shipRng)
% adds ship location from AIS data to a CSA plot based on the radial
% velocity from the AIS data. PLOTS AS A PATCH

a=axis;

hp=patch([min(shipRadVel) max(shipRadVel) max(shipRadVel) ...
          min(shipRadVel) min(shipRadVel) ],a([3 3 4 4 3]),'k');
set(hp,'EdgeColor',[.7 .7 .7],'FaceAlpha',0.2)


% h=patch([min(shipRadVel) max(shipRadVel) max(shipRadVel) min(shipRadVel) min(shipRadVel) ],a([3 3 4 4 3]),'c');
% set(h,'EdgeColor','c','FaceAlpha',0)
% plot(mean(shipRadVel),a(4)-(a(4)-a(3))/10,'k.')
% text(mean(shipRadVel).*1.1,a(4)-(a(4)-a(3))/10, ...
%     ['Ship Range: ' num2str(min(shipRng),3) ' to ' num2str(max(shipRng),3) 'km'])

end

function test_case


fn ='/m_files/test_data/cs_plot/CSS_cop1_06_05_30_2100.cs4';

CS = ReadCS(fn);


%m = 7 ; % CSQ_Rfg1_10_07_05_182228.cs
rcIdx = 3;
%fbinIdx = 202;

% Plot to veryfy .
figure
H = cs_plot(CS,rcIdx); hold on,


LS = line_style_groups;

h = cs_plot(CS,rcIdx,LS(2),H,160);

% try char inputs for line styles
h = cs_plot(CS,rcIdx,'r*',H,160);


keyboard

end


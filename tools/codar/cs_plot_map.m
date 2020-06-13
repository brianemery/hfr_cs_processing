function H = cs_plot_map(CS,peakIdx)
% CS PLOT MAP - 3D plot of CS data like SpectraPlotterMap
% S = cs_plot_map(CS)
% 
% INPUT
% CS        - cross spectra data structure
% peakIdx   - (optional) cell array of first order points
% 
% OUPUT
%
% EXAMPLE
%
% figure, H = cs_plot_map(CS)
% w
% hold on
% [freqs,~,~] = getDopplerVelocities(CS.Header);
% plot(H.axes(1),freqs(peakIdx{1}([1 end])), [ 1 1],'k.' )
% plot(H.axes(1),freqs(peakIdx{1}([1 end])), [ 1 1],'k*' )
% plot(H.axes(1),freqs(peakIdx{1}), ones(size(peakIdx{1})),'k.' )
% plot3(H.axes(1),freqs(peakIdx{1}), ones(size(peakIdx{1})), ones(size(peakIdx{1})),'k.' )
%
% This function desperately needs the ability to add the FOL's
% .. this is sort of close
% for j = 1:numel(peakIdx), plot3(H.axes(1),freqs(peakIdx{j}), j*ones(size(peakIdx{j})), ones(size(peakIdx{j})),'k.' ), end


% Copyright (C) 2014 Brian Emery
%
% Version 13-Aug-2014 12:40:01

% check for test case
if strcmp('--t',CS), test_case, return, end
   

% WHOI data has length(freqs) that is half what we'd expect due to trimming
% of the matricies, but the freqs array is in the CS file
if  isfield(CS,'freqs') && ~isempty(CS.freqs) %&& nargin < 2
    freqs = CS.freqs;
else
    [freqs,~,~] = getDopplerVelocities(CS.Header); 
end

% LERA specific contour range
if isfield(CS,'SpecHead')
   clvl = [-120 -50];
else
    clvl = [-165 -100]; %[-205 045])
end


% Get correct field name
fn = {'antenna3Self','a33'};
fn = char(fn(isfield(CS,fn)));


% get # range cells
rc = size(CS.(fn),2);

% convert to dbm
CS = cs_volts2dbm(CS);

% Define fields to plot
fn = {'antenna13CrossSp','antenna23CrossSp','antenna3Self', ...
      'a13','a23','a55'}; % was a33
fn = (fn(isfield(CS,fn)));


% and the labels
fnlab = fn; %{'A13','A23','A33'};


% get the first order region information
if nargin > 1
    
    ilim = peak_to_index(peakIdx, size(CS.freqs,1)/2 );
    
    [x,y] = alim_to_xy(freqs,1:rc,ilim);
    
else
    [x,y] = deal([]);
    
end



% MAKE FIGURES
figure 

hx = makesubplots(3,1,.05,.05);

for i = 1:3
    
    h = surf(hx(i),freqs,1:rc,CS.(fn{i}).');
    
    set(h,'EdgeColor','interp')
    
    view(2)
    
    caxis(hx(i),clvl) %[-205 045])
        
    % add FOL's?
    plot_fols(hx(i),x,y)
    
    hc = colorbar('peer',hx(i));
    
    ylabel(hx(i),'Range Cell')
    
    colorbar_label([fnlab{i} ' (dBm)'],hc)
       
    H.surf(i) = h;
    
end




set(hx(1),'XTickLabel','')
set(hx(2),'XTickLabel','')
xlabel(hx(3),'Doppler Frequency (Hz)')


ht = title(hx(1),[fileparts_name(CS.FileName)]);
set(ht,'Interpreter','none')


set(gcf,'Units','inches')
set(gcf,'Position',[12.7500    7.7222    9.4583    9.1250])

H.axes = hx;
H.title = ht;

end

function test_case
% TEST CASE
% 
% test case directory: /m_files/test_data/


load('/projects/hf_winds/data/cop1/csq/CSQ_cop1_20_02_01_000141.mat')

% FOL user parameters (see imageFOL.m help)
user_param = [40 100 10]; % [stdev maxcurr snrcutoff] <-- maybe use 10 for SNR? 

[~,Vrad,~] = getVelocities(CS.Header);

[~,~,dv] = getDopplerVelocities(CS.Header);

n = length(Vrad)/2;

% get approx indecies of Bragg lines
[~,iBragg(1)] = min( abs( Vrad(1:n)) ) ;
[~,iBragg(2)] = min( abs( Vrad(n:end)) ) ; 

% adjust last result
iBragg(2) = iBragg(2) + n - 1;

% velocity increment m/s
v_incr = dv/100; %mode(diff(Vrad))/100;

[FOreg, FOregi, Alims] = imageFOLs(CS.antenna3Self.',iBragg,v_incr,user_param);


% convert FOregi to peakIdx
peakIdx = convert_fol_ix(FOregi);




H = cs_plot_map(CS,peakIdx)


return 

% other possible tests

% csq file to look at
csqNm = '/m_files/test_data/getFirstOrder/CSQ_Rfg1_10_07_21_060515.cs';

% Read it in ...
CS = ReadCS(csqNm,1);

% Get frequencies
[CS.freqs,CS.Vrad] = getVelocities(CS.Header);


% Run the test
[peakIdx,Alims,absIdx,abstf] = getFirstOrder(CS,CS.Vrad);




end


% All code for adding first order lines
function ilim = peak_to_index(peakIdx,mx)
% Convert the cell array of indecies to something I can plot
%
% mx is the mid point index for splitting
%
% This function assumes that peakIdx has a cell for every range cell

ilim = NaN(length(peakIdx),4);

for i = 1:numel(peakIdx)
    
    px = peakIdx{i};
    
    % split left and right 
    left  = px( px < mx ); 
    right = px( px > mx );
    
    % create the outupts (which are the x indecies)
    if ~isempty(left)
        ilim(i,1:2) = [min(left)-1 max(left)+1];
    end
    if ~isempty(right)
        ilim(i,3:4) = [min(right)-1 max(right)+1];
    end

end


end

function [x,y] = alim_to_xy(freqs,rc,ilim)
% ALIM TO XY - 
% convert first order indecies to x and y for the plot
%
% this assumes that each row of ilim corresponds to a range cell

% Copyright (C) 2020 Brian M. Emery 
%
% The COVID-19 pandemic avoidance stay at home with the kids version

[x,y] = deal(NaN(size(ilim)));

% ... I guess use loops to avoid NaN?

for i = 1:size(ilim,1)
    for j = 1:size(ilim,2)
        if ~isnan(ilim(i,j))
            x(i,j) = freqs(ilim(i,j));
            y(i,j) = rc(i);
        end
    end
end




end

function plot_fols(h,x,y)

if ~isempty(x)
    plot(h,x,y,'k-')
end
   

end
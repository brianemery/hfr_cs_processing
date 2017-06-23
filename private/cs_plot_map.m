function H = cs_plot_map(CS)
% CS PLOT MAP - 3D plot of CS data like SpectraPlotterMap
% S = cs_plot_map(CS)
% 
% in
% 
% OUPUT
%
% EXAMPLE
%

% Copyright (C) 2014 Brian Emery
%
% Version 13-Aug-2014 12:40:01

% check for test case
%if strcmp('--t',CS), test_case(mfilename), return, end

if nargin < 2, [freqs,~,~] = getDopplerVelocities(CS.Header); end

% get # range cells
rc = size(CS.antenna3Self,2);

% convert to dbm
CS = cs_volts2dbm(CS);

% Define fields to plot
fn = {'antenna13CrossSp','antenna23CrossSp','antenna3Self'};

% and the labels
fnlab = {'A13','A23','A33'};



% MAKE FIGURES
figure 

hx = makesubplots(3,1,.05,.05);

for i = 1:3
    
    h = surf(hx(i),freqs,1:rc,CS.(fn{i}).');
    
    set(h,'EdgeColor','interp')
    
    view(2)
    
    caxis(hx(i),[-165 -100]) %[-205 045])
    
    % colorbar('peer',hx(i))
    
    ylabel(hx(i),[fnlab{i} ' (dBm)'])
    
    
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

S = new_mfile(in);


end


% EXAMPLE TESTING FROM get_file_list.m

function test_case2
% TEST CASE

% flist base value (CELL)
base = {'/m_files/test_data/move_file/src/test_file.txt'};

% check listing works
flist = get_file_list('/m_files/test_data/move_file/src','test*');

check_result(flist, base, 'Test 1')    
    

% check empties returned if no files found
flist = get_file_list('/m_files/test_data/move_file/dest/','test*');

check_result(flist, {}, 'Test 2')


% check multiple input directories
wd ={'/m_files/test_data/move_file/src', '/m_files/test_data/move_file/dest',};
flist = get_file_list(wd,'test*');

check_result(flist, base, 'Test 3')


% test a good sized list for time ...
flist = get_file_list({'/m_files/test_data/cs_filter/good_cop_case', ...
                       '/m_files/test_data/cs_filter/bad_cop_case'},'CSQ*');

                   
% load previous result and check
load('/m_files/test_data/get_file_list.mat','prev')

check_result(flist, prev, 'Test 4')


% check empty returned for no directory found
flist = get_file_list('/m_files/test_data/move_file/dest/none','test*');

check_result(flist, {}, 'Test 5')

keyboard
end

function check_result(flist, base, testStr)

if isequal(flist,base)
    disp([testStr ' ... ok'])
else
    disp([testStr ' ... NOT ok'])  , keyboard 
end


end

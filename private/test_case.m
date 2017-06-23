function test_case(tst)
% DOA TEST CASES - tests for the various DOA methods, in one place
%
% INPUT(S)
% fxn  - name of the function to test
% tst  - name of the test to run? (tonys, sim ...)
%
% would be great to just grab the name of the calling mfile .. I know i've
% done that before

% ADD THIS TO:
% mle.m
% mle_ap.m
% music.m

% F! too complicated, need to simplify this 

% POST DEVELOPMENT, THIS IS HOW THIS WILL WORK :::

% *** BROKEN! ***


% GET NAME OF CALLING FUNCTION
% from creation_info.m
[st,~]=dbstack('-completenames');  
[~,fun,~] = fileparts(st(end).file);



% create the default test?
if nargin == 0 , tst = 'tonys'; end


switch tst
   
    case 'tonys'
        
        tonys(fun)
    
    case 'sim'
        
        sim(fun)
        
end

end

% TEST CASES SPECIFIED HERE
function sim(fun)
% TEST CASE - test DOA methods using simulation data
% 
% 


% TO DO
% more bearings
% need a way to show which bearings the velocities are truely comming from
%   ie, something about binning the data ...


% ----------------------------------------------------------------------
% TEST USING CSSIM DATA
% Test data made with radar_simulation.m:
[CS,APM,S] = radar_simulation('--');

% load /projects/hf_doa_methods/test_data.mat APM CS simTh Vr

% Compute radial (ocean current) velocities
[~,CS.Vrad] = getVelocities(CS.Header);

% Get First order region
% peakIdx is a cell array with left and right indecies for each range cell.
[peakIdx,~,~,~] = getFirstOrder(CS,CS.Vrad);

% choose range cell
rdx = 1;
 

% PLOT TO CHECK GET FIRST ORDER

H = cs_plot(CS,1); hold on,
LS = line_style_groups;
h = cs_plot(CS,rdx,LS(1),H,peakIdx{1});



% % ADD TO PLOTS INCREMENTALLY
% % plot the input currents
% h = plot_currents(APM.BEAR,S.Vr); 


keyboard % specify function??

% RUN DOA ALGORITHMS
% should just output the bearing and velocity ...
R = doa_on_range_cell(CS,APM,peakIdx{rdx},rdx);



% PLOTS

% plot the input currents
h = plot_currents(APM.BEAR,S.Vr); hold on

% MUSIC single and dual
hu    = plot(R(1).Bear(:,2:3),R(1).RadVel(:,2:3),'go');
hu(1) = plot(R(1).Bear(:,1),R(1).RadVel(:,1),'g*');

% MLE single and dual
hx    = plot(R(2).Bear(:,2:3),R(2).RadVel(:,2:3),'ro');
hx(1) = plot(R(2).Bear(:,1),R(2).RadVel(:,1),'r*');

legend([hu(:)' hx(:)'],'MUSIC dual','MUSIC sngl','MLE dual','MLE sngl')

keyboard


% NEED HYPOTHESIS TESTING!!

return
% OLDER TO RECYCLE?

% add Bragg velocity, and bearings from MLE
h1 = plot(APM.BEAR(ix(:,1)),Vrad(fbin),'r*');          % ... single
h2 = plot(APM.BEAR(ix(:,2:3)),Vrad([fbin fbin]),'g*'); % ... dual

% add Bragg velocity, and bearings from MLE
h3 = plot(APM.BEAR(sdx),Vrad(fbin),'ro');          % ... single
h4 = plot(APM.BEAR(ddx),Vrad([fbin fbin]),'go'); % ... dual


keyboard




% ---- 
% CHECK PLOT 
% Make sure I know what is going on ...
% Should show dual bearing situation

% plot whole velocity profile
plot(simTh*180/pi,Vr), hold on 

% find everything contributing to this fbin
dv = mode(diff(Vrad))/2;
ix = find( Vr > (Vrad(fbin) - dv) & Vr < (Vrad(fbin) + dv) );

plot(simTh(ix)*180/pi,Vr(ix),'g*')

% ---- 


% TO DO : makes ure I'm using the right APM!!


end

function tonys(fun)
% TONYS TEST ... EXAMPLE FROM 
% Properties of HF RADAR Compact Antenna Arrays and Their Effect 
% on the MUSIC Algorithm by dePaolo and Terril
%
% DUAL! BEARING EXAMPLE FROM GETDOA.M
% 
% Makes a figure showing the DOA function vs bearing just like fig 9 in the
% De Paolo and Terrill Scripps report
%
% Lots of code from music.m
% 
% (APM,DOA,singleIdx,dualIdx)



% DATA FROM TONYS PAPER
[C,APM] = tonys_test_data;


% Display what's going on
disp(['TESTING ' upper(fun) '.M'])

% Create the function handle
fxn = str2func(fun);





% RUN THE DOA FUNCTIONS
[ix,TR] = mle_ap(APM,C,1);




keyboard


% ... when this is just a regular test ...
% ix = fxn(APM,C); 


% .... 

keyboard




% .... MORE TEST CODE ....

APM.BEAR(ix)

keyboard

% MUSIC SOLUTION
[DOA,singleIdx,dualIdx] = music(APM,C);

% instead of using getDOAs, use the fact that the known signal input
% bearing are 205 and 330 deg (from the reference)
dualIdx = [find(A.BEAR == 205) find(A.BEAR == 330)]; %APM?


% FIGURES
figure

subplot(211)
plot_doa(APM,DOA(:,1),singleIdx,'Single Bearing DOA function')
axis([0 360 -4 12])


subplot(212)
plot_doa(APM,DOA(:,2),dualIdx,'Dual Bearing DOA function')
axis([0 360 -5 30])


% DISPLAY CHECK NUMBERS
check_results(DOA(:,1),APM.BEAR,225)
check_results(DOA(:,2),APM.BEAR,205)
check_results(DOA(:,2),APM.BEAR,330)


keyboard
 	
% 
% 
% keyboard
% 



% ix = find(APM.BEAR == 330)
% DOA(ix,2)
% 10*log10(DOA(ix,2))
% ix = find(APM.BEAR == 205)
% 10*log10(DOA(ix,2))

end



% COMMON SUBFUNCTIONS HERE

function h = plot_currents(bear,Vr)
% PLOT CURRENTS

figure

h = plot(bear,Vr,'-b*'); hold on
% h2 = plot(bear,simVr1,'-b*'); 

xlabel('deg CWN'), ylabel('cm/s')
% legend([h1 h2],'Range Cell 1','Range Cell 31')

end

function plot_doa(APM,DOA,idx,titleStr)
%
% idx gives the location of the max of the DOA, single or dual angle
% solutions

plot(APM.BEAR,10*log10(DOA),'-b.')

hold on

plot(APM.BEAR(idx),10*log10(DOA(idx)),'g*')

xlabel('bearing (deg CWN)'),ylabel('10*log10(DOA)')

title(titleStr)

% add info to plot
for i = 1:length(idx)
text(APM.BEAR(idx(i))+10, 10*log10(DOA(idx(i))), ['(' num2str(APM.BEAR(idx(i))) ',' num2str(10*log10(DOA(idx(i)))) ')'])

end

end

function check_results(DOA,BEAR,brg)
% CHECK RESULTS
% 
% quickly check my resuls vs his

mine = num2str(real(10*log10(DOA(BEAR == brg))),4);

switch brg
    case 205
        his = '28.6';
        
    case 225
        his = '9.5';
        
    case 330
        his = '21.1';
end

disp(['DOA Metrics: his = ' his ' mine = ' mine ])


end
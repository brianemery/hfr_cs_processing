function  absix = mle_test(APM,C)
% MLE ALGORITHM - TESTING/DEV VERSION?
%
% C is the covariance matrix 

% TO DO
% What should this output?
% Test for single/dual? 
% 
% 
% probably efficiencies to be had pre-computing all the P matricies first
% and storing them?
%
% Other improvements based on refs that cite this method?
% More extensive study?
%
% Parse this out into functions, move music.m into hf_doa_methods?
%
% MLE/MUSIC outline:
% load/get CS
% get first order region (whole CS or range cell?)
% run doa (loop over range cells?)
% run single dual tests
% make the radial file 

% tic

% check for test case
if strcmp('--t',APM), test_case;  end

% SINGLE BEARING SEARCH
% Find the single bearing solution
% will probably need the value of the trace to compare with dual brg
% solution
[tr,ix] = find_mle(APM,C);



% TWO BEARING SEARCH
% now, combine this (one brg) soln with the other elements of the apm,
% compute and store the trace

% First compute the location of the 2nd source, assuming the first is at ix
% from above

[absix,maxtr] = deal(NaN(11,1));

% Initialize the (absolute) bearing index, and the values of the trace
absix(1) = ix;
maxtr(1) = tr;

% INIT THE ITERATIVE ALGORITHM
m = 2;
dtr = 1;

% for m = 2:10
while dtr > 10*eps

    % find the abs index of the APM to fix
    [tr,ix] = find_mle_dual(APM,C,ix);
    
    
    maxtr(m) = tr;
    absix(m) = ix;
     
    % check for convergence
    dtr = abs( real(maxtr(m)) -  real(maxtr(m-1)) );
    
    
    % Put a break in here?
    if m == 10,
        disp('MLE: 10 iterations reached, terminated early')
        break
    end
    
    % increment m
    m = m + 1;

        
end

% toc




% CREATE OUTPUTS
% output the indecies of the single and dual solutions

% strip off nan's?
absix = absix(~isnan(absix));

if length(absix > 3)
    absix = absix([1 end-1 end]);
else
    % ???
    keyboard
end


return
keyboard

% PLOT TO LOOK FOR CONVERGENCE
plot(1:numel(maxtr),real(maxtr),'-o')


% keyboard
% 
% [tr3,ix3,b3]= find_mle_dual(APM,C,ix2);
% 
keyboard
% 

% SINGLE VS DUAL (VS 3?) HYPOTHESIS TESTING



end


function [mxtr,ix] = find_mle(APM,C)
% FIND MLE - single bearing case

% Get number of bearings to loop over
nt = length(APM.BEAR);

tr = NaN(nt,1);

for n = 1:nt
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    % make it a column vector (M element x d signals), M = 3
    % A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
    
    % PRE CONSTRUCT THE APM MATRIX
    % If nargin == 2 we're looking for a 2 brg solution and need to modify the
    % APM search    
    A = [ APM.A13R(n)+ 1i.*APM.A13I(n)  APM.A23R(n) + 1i.*APM.A23I(n)  1+1i.*0 ].';

    
    
    % max(trace(P*C))
    % C is the covariance matrix (R in the //Introduction to DOA
    % Estimation// book)
    % A is the APM matrix with rows for each element (Mxd), d is signals? (pg 38)
    
    % Here A is A(theta), and P is a function of theta (thus loop over theta?)
    % do this in a subfunction to keep notation simple?
    % Note: ' is the Hermitian transpose, which is intended
    P = A*(inv(A'*A))*A'; % can I remove the inverse?
        
    tr(n) = trace(P*C);
        
    
end



% Ok ... what have I got?
[mxtr,ix] = max(tr); 

end

function [mxtr,ix] = find_mle_dual(APM,C,ii)
% FIND MLE - loop over APM bearings, finding the best fit to the signal

nt = length(APM.BEAR);

% TWO BEARING SEARCH
% create the APM excluding the ix bearing
B  = subsref_struct(APM,setdiff(1:nt,ii),nt,1);

% get the shortened index length
nt = length(B.BEAR);

tr = NaN(nt,1);

    
for n = 1:nt
    
    % Note that sqrt(sum(A.^2)) == sqrt(2) (Tony's section 4.2)
    % make it a column vector (M element x d signals), M = 3
    % A = [ APM.A13R(idx)+1i*APM.A13I(idx)  APM.A23R(idx)+1i*APM.A23I(idx)  1+1i*0 ].';
    
    % CONSTRUCT THE APM MATRIX
    % Construct the two signal version
    A = [ APM.A13R(ii)+ 1i.*APM.A13I(ii)  APM.A23R(ii) + 1i.*APM.A23I(ii)  1+1i.*0; ...
            B.A13R(n) + 1i.*B.A13I(n)       B.A23R(n)  + 1i.*B.A23I(n)     1+1i.*0; ].';
        
        
    
    % max(trace(P*C))
    % C is the covariance matrix (R in the //Introduction to DOA
    % Estimation// book)
    % A is the APM matrix with rows for each element (Mxd), d is signals? (pg 38)
    
    % Here A is A(theta), and P is a function of theta (thus loop over theta?)
    % do this in a subfunction to keep notation simple?
    % Note: ' is the Hermitian transpose, which is intended
    P = A*(inv(A'*A))*A'; % can I remove the inverse?
        
    tr(n) = trace(P*C);
        
    
end



% Ok ... what have I got?
[mxtr,iy] = max(tr); 


% % TEST PLOT
% figure
% plot(B.BEAR, real(tr)), hold on,
% plot(B.BEAR(iy),real(tr(iy)),'r*')


% ... so my two solutions so far are ix, and ii
% need to translate ix into the APM index

% TRANSLATE IX TO ABS INDEX
ix = find(APM.BEAR == B.BEAR(iy));


% TO DO
% look at A, P at the two bearings with really high tr and make sure the
% math is working right ...

end





function test_case
% TEST CASE 
% 
% 

% TO DO
% more bearings
% need a way to show which bearings the velocities are truely comming from
%   ie, something about binning the data ...


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







% -------------------------------

% ANOTHER TEST ... EXAMPLE FROM 
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


% Create the covariance matrix from de Paolo's example:
C=[ 0.2162          0.0303-0.0090i  0.3170-0.0063i; ...
    0.0303+0.0090i  0.0436         -0.0091+0.0213i; ...
    0.3170+0.0063i -0.0091-0.0213i  0.5416];

% Create the idealized pattern 
APM = make_ideal_pattern(225, 0:5:360);


ix = mle_test(APM,C); % !!! THIS WORKS
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
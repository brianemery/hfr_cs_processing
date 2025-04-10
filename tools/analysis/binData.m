function [avex,avey,N,stdev,stats]= binData(x,y,binWidth,binCenters,wts,isgauss)
% BIN DATA - Bin Average over y data using binWidth 
% [avex,avey,N,stdev,stats]=binData(x,y,binWidth,binCenters,weights,isgauss)
%
% Produces the average of x in each bin defined by y +/- YbinWidth/2. Thus,
% YbinWidth has same units as y. 
%
% Originally designed to bin over columns, ie y is a time variable, and
% columns of x are different times. 
%
% INPUTS:
% - x can be a matrix as long as the number of columns is equal to the 
%   length of y. (need to clarify this)
% - binWidth input must be scalar, or if specifying bin centers, must be
%   the same size as the binCenters vector. The binWidth and binCenters 
%   vectors can specify any bin width and center locations.
% - optional 5th input: weights for computing weighted average, same size
%   the x and y inputs
% - Gaussian window uses the binWidth/3 as the stdev of the Gaussian (see
%   subfunction).
%
% A final output, 'stats' is an optional structure containing the hi, low
% and median values in each bin.
%
%
% % EXAMPLE 1:
% given x = amplitude
%       y = bearing
% [aveX,aveY]=binData(x,y,5); produces 5 degree binned amplitude data
%
% [avex,avey,N,stdev]=binData(x,y,binWidth) can also be used to output the
% number of data points in each bin, and their standard deviation 
% 
% NaN's ignored throughout. Uses the find command, so no assumptions are
% made about data spacing (in time for example).
% 
% % EXAMPLE 2: 
% % define spatial separation bins as exp(x) with x=1.5:0.33:8.0. This gives
% % exponential bins:
% bc = [0 exp(1.5:0.33:8.0)/1000];
% bw=bc(3:end)-bc(1:end-2);
% bw=bc(3:2:end)-bc(1:2:end-2);
% bc = mean([bc(3:2:end); bc(1:2:end-2)]);
% [ASD_bin,VAR_bin,N,stdev]=binData(ASD(:),VAR(:),bw,bc);
%
% 
% For confidence intervals, use confidence_intervals.m
%
% SEE ALSO
% bin_like_hfr.m, bin_drifters_like_radials.m, binData.m,
% bin_data_struct.m, merge_rads_spatially.m, bin_radials.m
%
% FOR NOTES ABOUT RADIAL BINNING
% see /m_files/tools/codar/Contents.m 


% Copyright (C) 2009-2010 Brian Emery
% Nov 2009 - expanded
% Jul 2012 - added weights option, minor clean up

% TO DO
% - build a test case that checks the bin centers/widths works properly
%   (compare withe mv ave or something?)
% - see bin_apm_struct.m on a better way to organize the meta data struct
%   and unpack the stats substruct.


% --------------------------------------------------------- 
%  SORT OUT INPUTS
% --------------------------------------------------------- 


% test case
if strcmp(x,'--t')
    test_case, return
end

% eject if empty
if isempty(x), [avex,avey,N,stdev,stats] = deal([]); return, end

if size(x,2)== 1 || size(x,2) ~= length(y) 
    x = x'; 
end

% Check that the length of the binWidths is either 1, or the same size as
% the bin centers
if nargin < 4
    
    % check inputs
    if length(binWidth)>1
        disp('binWidth must be scalar')
        keyboard
    end
    
    if binWidth<1
        % using floor and ceil dont work with very small bins
        binCenters=(min(y)):binWidth:(max(y));
    else
        binCenters=floor(min(y)):binWidth:ceil(max(y));
    end
   
    % Keep this around for a while for bw compatibility
    % [avex,avey,N,stdev]=binData_old(x,y,binWidth);
    
elseif nargin == 4

    % check inputs
    if length(binWidth) == 1
        binWidth = binWidth * ones(size(binCenters));
        
    elseif length(binWidth)> 1 && ~isequal(length(binWidth),length(binCenters))
        disp('length(binWidth)~=length(binCenters)')
        keyboard
    end
end

% CHECK WEIGHTS INPUT
if nargin < 5
    wts =[];
    
elseif ~isempty(wts)
    
    % make sure it's a row vector
    wts = wts(:)';
    
%     % check size
%     if size(wts,1) ~= size(x,1)
%         % not sure what I was doing here ...
%         % wts = ones(size(x,1),1)*wts;
%         
%         % this is useful to try
%         wts = wts';
%     end
    
    
    if size(wts,2) ~= size(x,2)
        disp('wts and x must have equal columns')
        keyboard
    end
end

if nargin < 6
    % gaussian  weighting defaults to off
    isgauss = false;
end


[avex,avey,N,stdev,stats]=binData_new(x,y,binWidth,binCenters,wts,isgauss);

end

function [avex,bc,N,stdev,stats]=binData_new(x,y,binWidth,bc,wts,isgauss)
% BIN DATA NEW
% newer version of bin data for variable bin widths

% Brian Emery 4 Nov 2009

if length(binWidth)==1
    binWidth = binWidth*ones(size(bc));
end

% create empty matricies
[avex,stdev,N,hi,lo,med]=deal(NaN(size(x,1),length(bc)));
index ={};

for i=1:length(bc)
    
    % create the index of the data going into each bin
    idx = find(y >= (bc(i)-(binWidth(i)/2)) & y <= (bc(i)+(binWidth(i)/2)) );
    
    % compute the average and ancillary data
    [avex(:,i),stdev(:,i),hi(:,i),lo(:,i),med(:,i),N(:,i)] = stats_noNaN(x(:,idx));
    
    % compute the weighted mean if given weights
    if ~isempty(wts) && ~isgauss   
        avex(:,i) = mean_weighted(x(:,idx),wts(:,idx));
        
        % * ADD WEIGHTED STD HERE TOO *
        % Note that matlab's std appears to compute weighted STD but I cant
        % seem to get it to get the right answer .. which means I'm
        % treading on thin ice I think ...
        
    elseif ~isempty(wts) && isgauss && ~isempty(idx)
       
        % (rely on stats_noNaN to provide a NaN result if idx is empty)
        
        % Gaussian weigting applied to the weights in y direction ... this
        % is being coded for radial HFR data so some of the decisions may
        % not generalize. 
        avex(:,i) = mean_gaussian_weighted(y(idx),x(:,idx),wts(:,idx),bc(i),binWidth(i));
                
    end
    
    % save the indexing
    index{end+1} = idx;
end

% create the stats struct
vars ={'hi','lo','med','index'};
for i = 1:numel(vars)
    eval(['stats.(vars{i}) = ' vars{i} ';'])
end


% % CHECK PLOTTING CODE
% 
% LS.Color = 'r';
% LS.MarkerFaceColor = 'r';
% LS.LineStyle = 'none';
% LS.Marker = 'o';
% 
% % plot all the data
% plot(y,x,'-b.'), hold on
% 
% % plot data in the bin
% for i = 1:length(bc)
%     plot(y(index{i}),x(index{i}),'c.')
% end
% 
% % plot the binned data
% plot(bc,avex,LS)
% 
% keyboard

end

function avex = mean_gaussian_weighted(t,x,wts,bc,bw)
% MEAN GAUSSIAN WEIGHTED
%
% Sort of. This takes the weights (eg SNR) and then applies a gaussian to
% them based on the location of t relative to bc. Wts can be all ones
%
% Assumes that bw is the same units as t, eg, degrees (not indecies)
%
% Even if things are on the edge of the window and get really small weights
% this will still produce a similar result ... so yeah
% wts.*gwts(:)'
% 
% ans =
% 
%     0.0029    0.0004
% 
% avex = mean_weighted(x,wts.*gwts(:)')
% 
% avex =
% 
%    -5.2933
% 
% avex = mean_weighted(x,wts)
% 
% avex =
% 
%    -5.2933
% 

keyboard % I am not useing gausswin correctly

% also:
% "you expect that the Gaussian is essentially limited to the mean plus or
% minus 3 standard deviations, "
% stdev = (N-1)/(2*alpha);
% ... so take the window, divide be 3, use that for stdev and compute
% alpha: (in bw units, not N?)
%stdev = bw/3; alpha = 

% the 5 here needs to be adjustable. Note that it's not normalized
% if bw is not odd, the interp setup gets challenging
gw = gausswin(bw,5); % w = gausswin(N,alpha);

% get integer window
hwin = (bw-1)/2;
        
gwts = interp1( (bc-hwin):(bc+hwin) ,gw,t);

% row vectors here
avex = mean_weighted(x,wts.*gwts(:)');

end

function [avex,avey,N,stdev]=binData_old(x,y,binWidth)
% Old version for scalar binWidth and output y locations based on the
% inputs.

% Keep this around for a while for bw compatibility


if length(binWidth)>1, disp('binWidth must be scalar'), keyboard, end

if binWidth<1
    % using floor and ceil dont work with very small bins
    avey=(min(y)):binWidth:(max(y));
else
    avey=floor(min(y)):binWidth:ceil(max(y));
end

avex(1:length(avey))=NaN;
stdev(1:length(avey))=NaN;
N(1:length(avey))=NaN;


for i=1:length(avey)
    idx=find(y>(avey(i)-(binWidth/2)) & y<=(avey(i)+(binWidth/2)));
    %avex(i)=nanmean(x(idx));
    [avex(i),stdev(i),hi,lo,median]=stats_noNaN(x(idx));
    N(i)=length(idx);
end

% % Test plot:
% figure, plot(avey,avex,'o')
% hold on
% plot(Wplot,N,'.')

end

function test_case
% TEST CASE - fairly robust test case to check that this is working ...


% GAUSSIAN (OR OTHER WINDOW) % WEIGHTED BIN AVERAGE
%
% I want it to be super obvious if this is working so let's make an impulse
% function and average that
t = 1:180;
x = ones(size(t));
x(90) = 50;
x(5) = 50;

figure
plot(t,x,'-o'), hold on

%  bin centers
bc = 1:180;

% bin width is 5^o
bw = 5*ones(size(bc));


% compute means using weight code, with all ones as weights
% this checks execution
[avex,avet,N,stdev,stats]=binData(x,t,bw,bc,ones(1,size(x,2)));

plot(avet,avex,'-x')

% No compute means using weight code, with all ones as weights
% and a gaussian window

bw = 10*ones(size(bc));
[avex,avet,N,stdev,stats]=binData(x,t,bw,bc,ones(1,size(x,2)),1);

plot(avet,avex,'.',avet,avex,'-')

keyboard



% WEIGHTED BIN AVERAGE

% Create input data to bin average
t = 1:.01:180;
x = [3; 2; 1]*sin(2*pi*t/180);

figure
plot(t,x), hold on


% put the bin centers near zeros for the sin curve
bc = 1:10:180;


% bin width is 90^o
bw = 90*ones(size(bc));

% compute means using weight code, with all ones as weights
% this checks execution
[avex,avey,N,stdev,stats]=binData(x,t,bw,bc,ones(1,size(x,2)));

plot(avey,avex(1,:),'-bo')
plot(avey,avex(2,:),'-go')
plot(avey,avex(3,:),'-ro')

[avex,avey,N,stdev,stats]=binData(x,t,bw,bc);

plot(avey,avex(1,:),'-b*')
plot(avey,avex(2,:),'-g*')
plot(avey,avex(3,:),'-r*')

% if sum(sum(avex > 1e-10)) > 0
%     disp('test failed!'), keyboard
% else
%     disp('testing binData.m ... ok')
% end

keyboard



% EQUAL BIN WIDTHS 

% Create input data to bin average
t = 1:.01:1000;
x = [3; 2; 1]*sin(2*pi*t/180);

figure
plot(t,x), hold on


% put the bin centers near zeros for the sin curve
bc = 90:90:900;

% bin width is 90^o, such that the bin means should be near zero
bw = 90*ones(size(bc));

% compute means
[avex,avey,N,stdev,stats]=binData(x,t,bw,bc);

plot(avey,avex(1,:),'-bo')
plot(avey,avex(2,:),':g.')
plot(avey,avex(3,:),'rs')

if sum(sum(avex > 1e-10)) > 0
    disp('test failed!'), keyboard
else
    disp('testing binData.m ... ok')
end

keyboard

% % start with 90, should produce zeros (or numbers very close to zero
% bc = 45:90:900;
% bw = 90*ones(size(bc));
% [avex,avey,N,stdev,stats]=binData(x,t,bw,bc);
% 
% plot(avey,avex(1,:),'-bo')
% plot(avey,avex(2,:),':g.')
% plot(avey,avex(3,:),'rs')
% 
% if sum(sum(avex > 1e-10)) > 0, disp('test failed!'), keyboard, end
% 
% 
% keyboard

end
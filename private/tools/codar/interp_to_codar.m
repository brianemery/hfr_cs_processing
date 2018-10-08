function [u,v,er,TUV] = interp_to_codar(lon,lat,stime,TUV)
% INTERP TO CODAR - get hf derived velocity given postion and time arrays
% [u,v] = interp_to_codar(lon,lat,stime,TUV)
%
% Given arbitrary time and position arrays (lon,lat and serial time), this
% will output the codar derived velocities interpolated in time and space.
%
% assumes that the number of columns of lon lat equal the elements of stime
% (time increments column-wise). Rows of lon lat can be, for example,
% positions of different drifters.
% 
% from get_codar_vel_at_drifter_pos.m
%
% For drifters, input the mid-fix positions and times. 
%
% Optionally, include the TUV struct to use as a 4th input

% Copyright (C) 2009-2010 Brian Emery


% Brian Emery
% 23 Nov 2009

%% --------------------------------------------------------- 
%  SETTINGS
%----------------------------------------------------------

% if nargin < 4
%     drive = '/Volumes/CodarData/Data/totals/';
% end


%% --------------------------------------------------------- 
%  GET DATA AND INTERP
%----------------------------------------------------------

if nargin < 4
    % Get Codar Data
    TUV = get_hfr_totals(stime);
    
end

% get central grid location
loc = [mean_noNaN(lon(:))  mean_noNaN(lat(:))];

% convert gridd to km
[HFx,HFy] = lonlat2km(loc(1),loc(2),TUV.LonLat(:,1), TUV.LonLat(:,2));
    [x,y] = lonlat2km(loc(1),loc(2),lon,lat);

% % GRID CODAR DATA (3-D)
[X,Y,Ug,Vg]=grid_codar_data(HFx,HFy,TUV.U,TUV.V);
%[X,Y,Ug,Vg]=grid_codar_data(TUV.LonLat(:,1),TUV.LonLat(:,2),TUV.U,TUV.V);

% GRID Fit Error
try
    [Xe,Ye,Err] = my_griddata(HFx,HFy,TUV.ErrorEstimates(1).TotalErrors); 
catch
    Err = NaN(size(Ug));
    disp('* Error field not found *')
end

% expand X, Y, and make T, all same size as Ug
X = repmat(X,[1 1 length(TUV.TimeStamp)]);
Y = repmat(Y,[1 1 length(TUV.TimeStamp)]);

T = reshape(TUV.TimeStamp,1,1,length(TUV.TimeStamp));
T = repmat(T,[size(X,1) size(X,2) 1]);

% Create empty HF outputs
[u,v,er] = deal(x.*NaN);%

% INTERP TO DRIFTER TIMES/LOCATIONS
% loop over drifter rows, and time columns
for r = 1:size(x,1)
    for c = 1:size(x,2)
        try
            u(r,c) = interp3(X,Y,T,Ug,x(r,c),y(r,c),stime(c));
            v(r,c) = interp3(X,Y,T,Vg,x(r,c),y(r,c),stime(c));
           er(r,c) = interp3(X,Y,T,Err,x(r,c),y(r,c),stime(c));
        catch
            
            [u(r,c),v(r,c),er(r,c)] = deal(NaN); %
        end
    end
end

end

%% -----------------------------------------------------------------------
%  DATA LOADING
function TUV = get_hfr_totals(stime)
% GET HFR TOTALS

% subfunction to keep things clean

% allow for simple usages
if length(stime) ==1
    stime = [floor(stime*24)/24 ceil(stime*24)/24];
end

% Set angular difference criteria
ang_dif = [30 150];

% load from 'database'
[U,V,stime,gridd,Count,AngDiff,Err] = get_total_data(stime([1 end]));



% Apply Maskings
Err(AngDiff< ang_dif(1) | AngDiff>ang_dif(2)) = NaN;
  U(AngDiff< ang_dif(1) | AngDiff>ang_dif(2)) = NaN;
  V(AngDiff< ang_dif(1) | AngDiff>ang_dif(2)) = NaN;

% Land mask
kk = sbgrid_mask(gridd);
U(kk,:) = NaN;
V(kk,:) = NaN;
Err(kk,:) = NaN;


% Put in TUV Structure
TUV = standardize_totals(AngDiff,U,V,gridd,stime,Err);






return

% %% --------------------------------------------------------- 
% %  CUSTOM CODE FOR ALASKA TOTAL DATA
% %----------------------------------------------------------
% % load([CFG.drive 'tot_pws.mat']);
% % 
% % [U,V]=deal(NaN(size(TUV.U,1),length(CFG.runTimes)));
% % serialtime = CFG.runTimes;
% % 
% % % Find and output data only at intersecting times
% % [c,iOut,j]=intersect(round(TUV.TimeStamp*24)/24,round(CFG.runTimes*24)/24);
% % 
% % AngDiff=[];
% % U(:,j) = TUV.U(:,iOut);
% % V(:,j) = TUV.V(:,iOut);
% % gridd = TUV.LonLat;
% 
% % This is real easy to adapt to TUV format when future calls for it.
% load([CFG.drive 'tot_pws.mat']);
% 
% % Put on contiguous time base so all the matricies are identical 
% [stime,TUV.U]=timebase(TUV.TimeStamp,TUV.U,1/24);
% [TUV.TimeStamp,TUV.V]=timebase(TUV.TimeStamp,TUV.V,1/24);
% TUV.ProcessingSteps{end+1} = 'timebase.m';
% 
% % Find and output data only at intersecting times
% [c,iOut,j]=intersect(round(TUV.TimeStamp*24)/24,round(CFG.runTimes*24)/24);
% 
% %TUV.AngDiff=[];
% % Need to make timebase operate on substructures in TUV (such as the error
% % structure). could use subsref_struct_cols.m here 
% % TUV = subsrefTUV(TUV,':',iOut);
% TUV.U=TUV.U(:,iOut);
% TUV.V=TUV.V(:,iOut);
% TUV.TimeStamp=TUV.TimeStamp(iOut);
% 
% % % Check for and potentially fill gaps if they are very short (using interp)
% % I = findGaps(TUV.TimeStamp,sum(~isnan(TUV.U),1)',2/24);
% % 
% % keyboard
% % if I(4) < 3/24
% %     disp('possible to use interp here ...')
% %     keyboard
% % end
% 
% return

end

%% -----------------------------------------------------------------------
% LOCAL GRIDDING FUNCTION 
function [X,Y,Ug,Vg]=grid_codar_data(HFx,HFy,U,V)
% GRID CODAR DATA
% [X,Y,Ug,Vg]=grid_codar_data(HFx,HFy,U_,V_)
% Turn arrays in to matricies, gridded in km space, creating the grid 
% from the inputs HFx, HFy. This makes the code for calc pos drifter much simpler.

% Brian Emery 3 Nov 2009

% Code here only works if the grid is somewhat regular (as in the case of
% sbchan grid when converted to km)
[X,Y,Ug] = my_griddata(HFx,HFy,U);
[X,Y,Vg] = my_griddata(HFx,HFy,V);

% % PWS code. This puts the codar data on a regular grid
% [X,Y,Ug] = my_griddata(round(HFx),round(HFy),U);
% [X,Y,Vg] = my_griddata(round(HFx),round(HFy),V);


return
%% -----------------------------------------------------------------------
% NOTES
% -----------------------------------------------------------------------
% FOR the PWS data, I ran both the above code using my_griddata, and the
% code below. My version resulted in about 20% fewer data points, which I
% would have expected for griddata.m.

% Clean up the grid by Rounding to nearest 1/2 km
[X,Y]=meshgrid(unique(round(HFx)),unique(round(HFy))');
 
% [XI,YI,ZI] = GRIDDATA(X,Y,Z,XI,YI) also returns the XI and YI formed
% GRIDDATA does some interpolating, so some edge points are lost ..
for i = 1:size(U,2)
    Ug(:,:,i)=griddata(HFx,HFy,U(:,i),X,Y);
    Vg(:,:,i)=griddata(HFx,HFy,V(:,i),X,Y);
end

% I think I want to use interp for cases when the grid is not uniform, ie
% this code cleans up messy grids
% need some code to check no missing values ...
% keyboard

% % CHECK PLOTS
% hh=arrowplot(HFx, HFy,U(:,1),V(:,1),1/10);
% 
% hold on
% hh=arrowplot(X, Y,Ug,Vg,1/10,'r');
% plot(X,Y,'r.')
% plot(HFx, HFy,'k.')
% 
% [X2,Y2,Ug2,address] = my_griddata(HFx,HFy,U(:,1));
% [X2,Y2,Vg2,address] = my_griddata(HFx,HFy,V(:,1));
% hh=arrowplot(X2, Y2,Ug2,Vg2,1/10,'c');
% plot(X2, Y2,'c.')



end


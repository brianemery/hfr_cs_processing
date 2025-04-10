function [FOreg, FOregi, Alims]=imageFOLs(ant3_volt,iFBragg,v_incr,user_param,islera)
% IMAGE FOL - define first order region using image processing
% 
% [FOreg,FOregi]=imageFOLs(Data,iFBragg,v_incr,user_param,islera);
%
%  uses Marker-Controlled Watershed Segmentation (available in the Matlab
%   image processing tool box) to find and separate
%   areas of peak energy in monopole (ant3) that should be with in the
%   first order bragg region from troughs that should indicate the space
%   between the first and second order regions.
%   Returns all indices of spectral
%   points within the FOL region where snr>=5.
%
% calls:
%   -multiple functions from the MATLAB image and signal processing toolboxes.
%   -getpeaks.m
%
% INPUTS
%   data         --  cross spectra matrix (nRC x nFFT)
%   iFBragg      --  indices of the +/- Bragg frequencies
%   v_incr       --  velocity resolution of spectra (m/s)
%   user_param   --  [vel_scale max_vel snr_min];  
%   islera       --  boolean set to true for LERA data
%
% where:
%
%        vel_scale  -- velocity scale in cm/s used to set N, the core smoothing
%                        lengthscale used throughout imageFOLs.  This can be defined
%                        either visually as 1/2 the spectral width of most
%                        of the Bragg energy or mathematically as the temporal std
%                        dev of the observed velocities.
%
%                        Typical values:
%                        velD_change=20; for 25 MHz
%                        velD_change=40; for 13 or 4-5 MHz
%
%        max_vel -- in cm/s: the absolute maximum current speeds ever encountered
%                       suggested value= 200.
%
%        snr_min -- minimum signal to noise ratio allowed in outgoing FOreg
%                       suggested value= 5.

%
%   !!! additional constants are set internally but do not normally vary with
%             site/conditions  !!!
%
%
% OUTPUTS
%          FOReg --  returning matrix (same size as data.a3) marking
%                       FO region spectra with snr>=5 (as 1, else 0)
%          FORegi -- range/dop velocity indices of all points marked as FO
%          Alims  -- a range cell x 4 matrix that gives the first order
%                       limits, as defined by COS.
%
%
% EXAMPLE
% 
%
%
%   Notes:
%   -Includes tests as applied in companion paper submitted to J. Atmos. and
%   Oceanic Techno. titled: 
%       "Improved Detection of First Order Region for Direction Finding
%        HF Radars using Image Processing Techniques"  by Anthony Kirincich 
%
%
%    Anthony Kirincich
%    WHOI-PO
%    akirincich@whoi.edu
%
% Version: _v07122016
%
% 2017 04 10 Modification by Brian Emery
%            Minor refactoring 
% 2019 09 26 Line 245 added if thens to handle empty values
% 2019 10 25 Added code/functionality from lera_imageFOLs_v5.m so this can now handle
%            both seasonde and lera. Created a test case for this in the
%            function imageFOLs_dev_test.m 

if nargin < 5, islera = false; end

% Set plotting switches
goplot = [0 0 ];


%%%%%%%%%%%% User-set parameters for processing %%%%%%%%%%%%
vel_scale=user_param(1);
max_vel=user_param(2);
snr_min=user_param(3);



%%%%%%%%%%%% set sFOL standard parameters for processing %%%%%%%%
% generally do not vary with site   
min_dB=-170;  %in dB, threshold lowest power.
seg_thres=4;   %min # of resulting MCWS segments that must be found in each half
% if lower # is found, Alim output is based on snr
% exceedances rather than segments alone  (see below).

%%%%% set the minium value of DN that is allowed.  %%%%%%%%%
%%% an example of a hard threshold
%min_DN=5;  
%%% an example of a dynamic threshold %%%%%%%
min_DN=ceil(.25*vel_scale/(v_incr*100)); %1/4 of the minimum smoothing lengthscale

%coefficients to filter the spectra for high noise levels that exist within the Bragg region itself
%  due to external RF noise that would normally pollute the FOL results.  Coefficients are set to
Bragg_noise_flag=[7 .444];
%   (1)  set a minimum value for the 'rolloff' of Bragg energy with range, estimated as the difference in the
%        max and min range dependent 'mean'  Bragg energy.  (see the rhs panel in Figure 2 for  details)
%   (2) set a minium value for the ratio between the rolloff of Bragg energy
%         for the left and right spectral halfs.   If the value is less than
%         the flag, the rolloff in the 'other side' is significantly greater,
%          suggesting a noise problem exists in this half with range.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%% Prep cross-spectra %%%%%%%%%%%%%%%%%%%%%%%%%
%%%    -estimating and subtracting range-dependent background energy
%%%    -normalizing by the range-dependent max energy

if ~islera
    
    % start by estimating the gain
    gain3 = 10*log10(abs(ant3_volt)) + (-40. + 5.8); 
        
    % set minimium value for gain3
    i=find(gain3<min_dB); gain3(i)=min_dB;
    
else
    % code here derived from lera_imageFOLs_v5 but for 2d matricies
    gain3 = 20*log10(abs(ant3_volt));
    
end

[n,m] = size(gain3);



%%%%% filter cross-spectra to smooth gain3  %%%%%%%
%initially, use a lengthscale based on the spectra length of vel_scale
%velocity change using v_incr
N=round(vel_scale/(v_incr*100));   %size of typical vel change
if (N/2)~=round(N/2) %N is odd reduce by 1
    N=N-1;
end

% set the maximum lengthscale
N_max=round(max_vel/(v_incr*100));   %size of the max vel change
if (N_max/2)~=round(N_max/2) %N is odd reduce by 1
    N_max=N_max-1;
end



%at this point N does not significantly effect the calculation, but does
%provide a measure of what the size (width) of bragg region should look like at
%this frequency.

%remove center frequencies... in case there is significant energy here and
%also provides a inner boundary for the MCWS calculation
ci=round(m/2-N):round(m/2+N);  %indices +/-N from the center freq.
gain3(:,ci)=nan;

%%%
%set up left and right hand side indices
li=1:round(m/2);  ri=round(m/2):m;

%get max, mean, and mean background for each range cell for each half
lmm=[nanmax(gain3(:,li)')'  nanmean(gain3(:,li)')' nanmean(gain3(:,li(1:N))')'];
rmm=[nanmax(gain3(:,ri)')' nanmean(gain3(:,ri)')' nanmean(gain3(:,ri(end-N+1:end))')'];

%%% filter result with range, using n/N as a lengthscale
f_l=ceil(n/N);
if (f_l/2)==round(f_l/2); %N is even, make odd
    f_l=f_l+1;
end
filt_shape=(ones(f_l,1)./f_l);   %create filter shape using square filter
%%%  filter for each side separately
for ii = 1:2
    if ii==1; g=lmm;    elseif ii==2; g=rmm;    end
    %make the filter form correct by padding ends with 1,end index pre conv.
    gg=[ones(floor(f_l/2),1)*g(1,:); g; ones(floor(f_l/2),1)*g(end,:)];
    h=[];
    for i=1:length(gg(1,:));
        h(:,i)=conv(gg(:,i),filt_shape,'valid');  %filter
    end
    %resave out to lmml or rmml
    if ii==1; lmml=h;    elseif ii==2; rmml=h;    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ID high noise within Bragg region by comparing the roll off of energy with range
Bragg_en=[mean(gain3(:,iFBragg(1)-N:iFBragg(1)+N),2) mean(gain3(:,iFBragg(2)-N:iFBragg(2)+N),2)];
g=Bragg_en;
gg=[ones(floor(f_l/2),1)*g(1,:); g; ones(floor(f_l/2),1)*g(end,:)];  %pad ends with 1,end index pre conv.
lBragg_en=[conv(gg(:,1),filt_shape,'valid') conv(gg(:,2),filt_shape,'valid') ];
Bragg_en_range=[max(lBragg_en,[],1)-min(lBragg_en,[],1)];   %the range of mean Bragg energy (BE update to apply to cols)
Bragg_en_ratio=[Bragg_en_range(1)./Bragg_en_range(2) Bragg_en_range(2)./Bragg_en_range(1) ];

%%% plot result of this RC-dependent function, if goplot(2)==1
if goplot(2)==1
    figure(8); clf;
    plot(lmm); hold on; grid on; plot(lmml(:,1:3),'k--');
    plot(rmm); hold on; grid on; plot(rmml(:,1:3),'k--');
    %plot(lmml(:,4),'g'); plot(rmml(:,4),'g')
end

%%% make h1: gain3-background energy state
h1=0.*ones(n,m);
h1(:,li)=gain3(:,li)-[lmml(:,2)*ones(1,length(li))  ];
h1(:,ri)=gain3(:,ri)-[rmml(:,2)*ones(1,length(ri))  ];

%%%% make h2: by normalizing h1 by max power for each half
h2=h1;
h2(:,li)=h2(:,li)./[(lmml(:,1)-lmml(:,2)) * ones(1,length(li))  ];
h2(:,ri)=h2(:,ri)./[ (rmml(:,1)-rmml(:,2)) *ones(1,length(ri))  ];

%remake center energy to 0s
h2(:,ci)=0;
i=find(h2(:)<0); h2(i)=0;
i=find(h2(:)>1.5); h2(i)=0;   % a weird fix for no power anywhere in the range cell, max = min and h2 blows up.
i=find(isnan(h2)==1); h2(i)=0;


%%%%%%%%%%%% Recast the lengthscale based on the 2nd order energy present %%%%%%%%%%
%%%
%%% Here, N_factor is defined, separately for each half of the spectra
%%%    as the difference between the mean 2nd order energy and the mean noise floor
%%%
%%%  N_factor is larger for larger waves, or 2nd order energies

%%%% old way, where 2nd order region location is static
% N_factor=[  mean(mean(gain3(:,[ (1:iFBragg(1)-3*N) (iFBragg(1)+3*N:ci(1)-1)] ))) - lmml(n,3) ;
%      mean(mean(gain3(:,[ (ci(end)+1:iFBragg(2)-3*N) (iFBragg(2)+3*N:m)] ))) - rmml(n,3) ]';

% ---- BE: this code deals with both halves of the spectrum ----

%dynamically define were the range-mean 2nd order region is, by looking for
%  the mean inner edge of the trough between 2nd and first order.
mm=mean(gain3);
%[pks,dzdt] = pksfinder(mm,ceil(2*nanstd(mm)));   %use the std dev as the threshold, helps for noise spectra
[pks,dzdt] = pksfinder(mm,ceil(1*nanstd(mm)));   %use the std dev as the threshold, helps for noise spectra
if goplot(2)==1
    figure(10); clf;
    plot(mm); hold on; plot(pks,mm(pks),'r+'); title('Range averaged spectra, to look for start of 2nd order')
    
end

% maxs are the tips of the bragg peaks
% old way
% maxs=[find(mm==max(mm(1:ci(1)))) find(mm==max(mm(ci(end):m)))];
% new that looks for the peak within bragg+/-N (from lera_imageFOL_v5.m, by BE)
maxs=[find(mm==max(mm(iFBragg(1)-N:iFBragg(1)+N))) find(mm==max(mm(iFBragg(2)-N:iFBragg(2)+N)))];


% find the peak that is to the left and right of each bragg peak, provided
% it is located in the correct quarter of the spectrum, correct if not.
pks2=nan.*ones(1,4); 

% --- BE: added two if then statements here to deal with empty i's. 
i = find(pks<maxs(1));               if isempty(i), pks2(1)=1;         else, pks2(1)=pks(i(end)); end
i = find(pks>maxs(1) & pks<ci(1));   if isempty(i), pks2(2)=ci(end)+1; else, pks2(2)=pks(i(1));   end
i = find(pks<maxs(2) & pks>ci(end)); if isempty(i), pks2(3)=ci(end)+1; else, pks2(3)=pks(i(end)); end
i = find(pks>maxs(2));               if isempty(i), pks2(4)=pks(end);  else,  pks2(4)=pks(i(1));  end  

% define the distance, in indices from the bragg peaks to each of the inner
% edges of the defined second order regions.
dpksbragg=pks2-[iFBragg(1)*ones(1,2) iFBragg(2)*ones(1,2)];

% define the start (outer) and ends (inner points) of the 2nd order spectrum
% around the Bragg peaks
i=find(abs(pks2-[iFBragg(1)*ones(1,2) iFBragg(2)*ones(1,2)])<3*N);
starts=[1 ci(1)-1 ci(end)+1 m];
ends=[iFBragg(1)-3*N iFBragg(1)+3*N iFBragg(2)-3*N iFBragg(2)+3*N];
starts(i)=ends(i);
ends(i)=pks2(i);

% okay, now define what the noise factor is as the difference between the
% 2nd order energy and the noise floor
N_factor=[  mean(mean(gain3(:,[(starts(1):ends(1))  (ends(2):starts(2))] ))) - lmml(n,3) ;
            mean(mean(gain3(:,[ (starts(3):ends(3)) (ends(4):starts(4))] ))) - rmml(n,3) ]';

%N_factor_scale=4./sqrt([min(abs(dpksbragg(1:2))) min(abs(dpksbragg(3:4)))]);   %old arbitrary method
%N_factor_scale=[mean(abs(dpksbragg(1:2))) mean(abs(dpksbragg(3:4)))];
N_factor_scale=[min(abs(dpksbragg(1:2))) min(abs(dpksbragg(3:4)))];

%use a banket adjustment to scale up N_factor if the 2nd order energy is closer than 3*N
i=find(N_factor_scale<N*3);  N_factor_scale(i)=2;
i=find(N_factor_scale>=N*3); N_factor_scale(i)=1;

%%% adjust velocity length scale (N) by N_factor to
%%%    narrow filter in higher wave/2nd order energy conditions.
%%%   The additional of N_factor_scale acts to further decrease N if the
%%%   2nd order energy is close to the Bragg as a big swell would create
DN=round([N - N_factor.*(N_factor_scale)]);

%%%   Note that at this point DN can really small or be negative, if the conditions are
%%%   crappy. Near 0 or negative values are not really viable to continue
%%%   with, and rather then install an additional arbitrary non-linear
%%%   scalar to prevent this, we'll use a lower limit to ensure the
%%%   smoothing length scale has SOME positive value able to filter small
%%%   scale noise in the MCWS process.
%%%
% insert min value of DN if appropriate
i=find(DN<=min_DN); DN(i)=min_DN;

%show what values you are moving forward with
disp(sprintf('ImageFOLs (DN N_factor): %2.0f %2.0f %2.2f %2.2f\n',DN,N_factor.*(N_factor_scale)))

%%% start plot of final result, if goplot(1)==1
if goplot(1)==1 
    figure(7); clf;
    subplot(211); pcolor(gain3); shading flat; colorbar; caxis([-160 -80]); title('Raw Ant3 Spectral Power with segments (red) and sFOLs (white) shown')
    subplot(212); pcolor(h2); shading flat; colorbar; caxis([0 1]); title('Normalized Ant3 Spectral Power with segments (red) and sFOLs (white) shown')
    num_plots=2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
%%%%%%%%%%%% Start MCWS processing %%%%%%%%%%

%clear outgoing variables
FOregi=[]; FOreg=[]; DL=[]; Alims=[]; Isave=[];
dn_save=[nan nan];

%%% loop over each half of the spectra and compute segments separately.
for ii=1:2
    if ii==1; si=li;
    elseif ii==2;  si=ri;      
    end
    H=h2(:,si);    %%%%% using h2!
    skipthisone=0; %reset skipping flag
    
    if Bragg_en_range(ii)>Bragg_noise_flag(1) & Bragg_en_ratio(ii)>Bragg_noise_flag(2)   %new flag for Bragg noise
             
        %more waves need a tighter disk, no waves, strong currents need a wide disk
        %allow for the lengthscale to decrease if it is too big (the smearing kills all the energy)
        igood=[]; count=0;
        while  isempty(igood)==1
            I=H; %transfer the data to I 
            count=count+1;
            dn=round(DN(ii)/count); %allow the lengthscale to reduce further if necessary
            
            %saturate the top percentile of the rc energy...
            %scale up p with rc, using i/dn to saturate a higher fraction at
            %high rcs
            for i=1:n;
                p=100-(round(2*dn+i/(2*dn)));
                if p<10; p=10; end
                p2=prctile(H(i,:),p);
                if p2<=0; p2=0.01; end;
                I(i,:)=(H(i,:)./p2); i=find(I(:)>1); I(i)=1;
            end
            
            %Create flat maxima inside each "object" using opening-by-reconstruction
            %and closing-by-reconstruction" to clean" up the image.
            %set the size of the dilating disk, to decrease if too big initially
            se = strel('disk',dn);
            % erode, reconstruct, dilate, and reconstruct the image
            Ie = imerode(I, se);
            Iobr = imreconstruct(Ie, I);
            Iobrd = imdilate(Iobr, se);
            Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
            Iobrcbr = imcomplement(Iobrcbr);       
            
            igood=find(Iobrcbr>0);
            if count>4;
                igood=1;
                skipthisone=1;
            end
        end %while statement
        %
        if skipthisone==0
            
            %Use the Sobel edge masks, imfilter, and some simple arithmetic
            %to compute the gradient magnitude. The gradient is high at the
            %borders of the objects and low (mostly) inside the objects.
            hy = fspecial('sobel');
            hx = hy';
            Iy = imfilter(double(Iobrcbr), hy, 'replicate');
            Ix = imfilter(double(Iobrcbr), hx, 'replicate');
            gradmag = sqrt(Ix.^2 + Iy.^2);
            
            %normalize gradmag and scale by the magnitude of smeared I matrix
            gm=(gradmag./max(gradmag(:))).*max(Iobrcbr(:))./1;
            Iobrcbr=Iobrcbr-gm;  %difference gm from Io to sharpen edges of major features
            
            Isave(:,si)=Iobrcbr;
            
            %find the regional maxima
            fgm = imregionalmax(Iobrcbr);
            %get rid of all maxima with # pixels smaller than the lengthscale
            bw2=bwareaopen(fgm,floor(N));
            %  prep and perform the segmentation
            D1 = bwdist(bw2);
            DL(:,si) =  watershed(D1)+((ii-1)*100);  %pack result in DL matrix
            i=find(DL==((ii-1)*100)); DL(i)=0;
            
            if goplot(2)==1;
                figure(8+ii-1); clf
                subplot(511); imshow(I); set(gca,'ydir','normal')
                if ii==1; title('lhs'); else title('rhs'); end
                subplot(512); imshow(Iobrcbr); set(gca,'ydir','normal')
                     title('Opening-closing by reconstruction (Iobrcbr)')
                subplot(513);  imshow(gm); set(gca,'ydir','normal')          
                subplot(514); imshow(Iobrcbr-gm); set(gca,'ydir','normal')        
                subplot(515);     I2 = I;    I2(bw2) = 255;
                     imshow(I2); set(gca,'ydir','normal') 
                     title('Regional maxima superimposed on original image (I2)')
            end
            %    Save dn
            dn_save(ii)=dn;

        else
            disp(['skipping the ' num2str(ii) ' half for lack of spectral energy'])
            DL(:,si)=nan.*H;
        end
        
    else
        disp(['skipping the ' num2str(ii) ' half for noise in the Bragg region'])
        DL(:,si)=nan.*H;
    end  % if skipping for Bragg noise
    
end  % ii loop over left and right halves
%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % post processing %%%%%%

%isolate the segments that intersect the bragg peak 
a=ceil(N/4);
useg=unique(DL(:,[iFBragg(1)-a:iFBragg(1)+a iFBragg(2)-a:iFBragg(2)+a])); nums=[];
i=find(useg>0); useg=useg(i);
is=[];
for i=1:length(useg)
    is=[is; find(DL(:)==useg(i))];
end
%break out into arrays of range cell and doppler vel indices
[x,y]=meshgrid(1:length(DL(1,:)),(1:length(DL(:,1)))');  %get grids of indices for this half
FOregi_all=[y(is) x(is)];
%get outline of FOL region and add in mag
FOreg_all=0*DL; %
FOreg_all(is)=DL(is); FOreg_all(is)=1;

%ignore results inside FOreg with snr<snr_min
i=find(h1(is)>=snr_min);
% FOregi=FOregi_all(i,:);
%get outline of FOL region and add in mag
FOreg=0*DL; FOreg(is(i))=1;


%count number of segments in each half, if lower than a threshold,
%use FOreg instead of FOreg_all to bound Alims
seg_no=[length(unique(DL(:,1:ci(1)))) length(unique(DL(:,ci(end):end)))];

%test whether there is only one segment within a spectral half,
%if so replace FOreg_all with FOreg for that half
if  seg_no(1)<seg_thres;
    FOreg_all(:,1:ci(1))=FOreg(:,1:ci(1));
end
if seg_no(2)<seg_thres;
    FOreg_all(:,ci(end):end)=FOreg(:,ci(end):end);
end

%whatever is used for FOreg_all, make sure the edges are zeros.
FOreg_all(:,[1:iFBragg(1)-N_max iFBragg(1)+N_max:iFBragg(2)-N_max iFBragg(2)+N_max:end])=0;

% Code from lera_imageFOLs_v5.m :
% just the points with SNR>5
FOreg(:,[1:iFBragg(1)-N_max iFBragg(1)+N_max:iFBragg(2)-N_max iFBragg(2)+N_max:end])=0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate output into something equvalent to COS's Alims
%define the outer edge of the bragg region as the first spectral point
% inside of the Bragg region segment. Bounded by iFBragg+/-N_max

dd=diff(FOreg_all')';
lhs=iFBragg(1).*ones(n,2); rhs=iFBragg(2).*ones(n,2);


% BE comment: this is where SS and LERA code deviate significantly, which
% is why I put these into subfunctions
% 
% This code loops over range cells to use information about the 2nd order
% that is in the close-by range cells

if islera
    
    [lhs,rhs] = set_outer_limits_lera(dd,iFBragg,lhs,rhs,n,m,N_max,FOreg);
    keyboard
    
else

    [lhs,rhs] = set_outer_limits_ss(dd,iFBragg,lhs,rhs,n,m,N_max,N);

end

%redefine outputs....
%redefine lhs and rhs as Alims
Alims=[lhs rhs];

%redefine FOreg_all and FOregi,  FOreg
FOreg_all=0.*DL;
for ii=1:n
    FOreg_all(ii,[Alims(ii,1):Alims(ii,2) Alims(ii,3):Alims(ii,4)])=1;
end

is=find(FOreg_all==1);
FOregi_all=[y(is) x(is)];

%ignore results inside FOreg with snr<snr_min
i=find(h1(is)>=snr_min);
FOregi=FOregi_all(i,:);

%get outline of FOL region and add in mag
FOreg=0*DL; FOreg(is(i))=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add to existing plot
if goplot(1)==1
    figure(7);
    for ii=1:num_plots
        subplot(num_plots,1,ii); hold on;
        [cs,ho]=contour(DL,[0 0],'r');
        [cs,ho]=contour(FOreg,[1 1],'w');
        plot(Alims,(1:n)'*[1 1 1 1],'k','linewidth',2)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % record this action to the HEAD
% HEAD.ProcessingSteps{end+1}=mfilename;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return

end

function [lhs,rhs] = set_outer_limits_ss(dd,iFBragg,lhs,rhs,n,m,N_max,N)
    
for ii=1:n
    
    % special treatment for first RC
    if ii==1
        i = [0 find(dd(1,:)~=0) m];
        
        ifl = find(i<iFBragg(1));       
        ifr = find(i>iFBragg(1));
        
        if abs(i(ifl(end))-iFBragg(1))<N_max; lhs(ii,1)=i(ifl(end)); end
        
        if abs(i(ifr(1))-iFBragg(1))<N_max; lhs(ii,2)=i(ifr(1)); end  
        
        ifl=find(i<iFBragg(2));     
        ifr=find(i>iFBragg(2));
        
        if abs(i(ifl(end))-iFBragg(2))<N_max; rhs(ii,1)=i(ifl(end)); end
        
        if abs(i(ifr(1))-iFBragg(2))<N_max; rhs(ii,2)=i(ifr(1)); end
        
    else 
        
        %handle left side of peak first.
        i=find(dd(ii,:)==1);
        
        if isempty(i)==0
            
            [s,is]=sort(abs(i-lhs(ii-1,1)));
            lhs(ii,1)=i(is(1));
            
            [s,is]=sort(abs(i-rhs(ii-1,1)));
            rhs(ii,1)=i(is(1));
        end
        
        %handle right side of peak second.
        i=find(dd(ii,:)==-1);
        
        if isempty(i)==0
            
            [s,is]=sort(abs(i-lhs(ii-1,2)));
            lhs(ii,2)=i(is(1))+1;
            
            [s,is]=sort(abs(i-rhs(ii-1,2)));
            rhs(ii,2)=i(is(1))+1;
            
        end        
        
        % handle cases with no real data
        if diff(lhs(ii,:))<N/4   
            lhs(ii,:)=iFBragg(1);
        end
        
        if diff(rhs(ii,:))<N/4   
            rhs(ii,:)=iFBragg(2);
        end
        
    end 
end



end

function [lhs,rhs] = set_outer_limits_lera(dd,iFBragg,lhs,rhs,n,m,N_max,FOreg)
%
% Code from lera_imageFOLs_v5
%
% ... many redundancies in this that should be simplified ...

for ii=4:n  %start this at RC 4 to be removed from wave returns in first RCs
    %%% if in the first RC, be careful to set the limits as close to the
    %%% Bragg line as possible to limit thes potential for 2nd-order swell
    %%% inclusion at low range cells,  i.e. a ridge line exists along the 
    %%% outer edge of the bragg, but curves back down to low RCs to include
    %%% 2nd order  (mostly a 25 MHZ problem.)
    if ii==3
        i=[0 find(dd(4,:)~=0) m];  %start this at rc3 not 1
        %For the left side 
        ifl=find(i<iFBragg(1));       ifr=find(i>iFBragg(1));    %find the ridges to the left and right of the bragg line
         fs=[(i(ifl(end))-iFBragg(1)) (i(ifr(1))-iFBragg(1))];
        if fs(1)>=-N_max & fs(1)<=0; lhs(ii,1)=i(ifl(end)); end   % place the left limit at the first ridge to the left of the bragg line that is less then N_max
        if fs(2)<=N_max & fs(2)>=0; lhs(ii,2)=i(ifr(1)); end       % place the right limit at the first ridge to the right of the bragg line that is less then N_max

        % for the right side
        ifl=find(i<iFBragg(2));     ifr=find(i>iFBragg(2));
         fs=[(i(ifl(end))-iFBragg(2)) (i(ifr(1))-iFBragg(2))];
        if fs(1)>=-N_max & fs(1)<=0; rhs(ii,1)=i(ifl(end)); end   % place the right limit at the first ridge to the right of the bragg line that is less then N_max
        if fs(2)<=N_max & fs(2)>=0; rhs(ii,2)=i(ifr(1)); end     % place the right limit at the first ridge to the right of the bragg line that is less then N_max
   %%% for the rest of the RC's 
    else  %ii>1
        %find the beginning and end of each segment at this range cell
        istart=find(dd(ii,:)==1); iend=find(dd(ii,:)==-1);
        if length(istart)~=length(iend) ; disp('segment error'); sdfasfasfad; return; end 
        %eliminate those segments that have no energy in the bragg region
        for i=1:length(istart)
            if isempty(intersect(istart(i):iend(i),find(FOreg(ii,:)==1)))==1;
                istart(i)=nan; iend(i)=nan;
            end
        end        
    istart=istart(~isnan(istart)); iend=iend(~isnan(iend));   %condense

      %%% now proceed to find the segments that overlap best        
        %handle left side of each of the peaks first. 
        if isempty(istart)==0
            [s,is]=sort(abs(istart-lhs(ii-1,1)));   %for the lhs
             if abs(istart(is(1))-iFBragg(1))<=N_max; lhs(ii,1)=istart(is(1)); end  %pick the one that is closest to the previous and within Nmax of the Bragg line
            %fix this estimate if the whole FOL is on one side of the Bragg
            %and there is a second segment that wraps thru bragg that I'm missing
             if istart(is(1))-iFBragg(1) > 0;  %it should be odd to have the starting value be on the wrong side of the bragg line
                if is(1)>1 %there is a start to the left of this one
                    if istart(is(1)-1)-iFBragg(1)<0 %the istart to the left of the shortest path is on the correct side of the Bragg line
                        lhs(ii,1)=istart(is(1)-1);   %use this one instead
                    end
                end
             end
             
            [s,is]=sort(abs(istart-rhs(ii-1,1)));   %for the rhs
             if abs(istart(is(1))-iFBragg(2))<=N_max; rhs(ii,1)=istart(is(1)); end
             %fix this estimate if the whole FOL is on one side of the Bragg
            %and there is a second segment that wraps thru bragg that I'm missing             
             if istart(is(1))-iFBragg(2) > 0  %it should be odd to have the starting value be on the wrong side of the bragg line
                if is(1)>1 %there is a start to the left of this one
                    if istart(is(1)-1)-iFBragg(2)<0 & abs(istart(is(1)-1)-iFBragg(2))<N_max %the istart to the left of the shortest path is on the correct side of the Bragg line
                        rhs(ii,1)=istart(is(1)-1);   %use this one instead
                    end
                end
             end
        end
        
        %handle right side of peak second.
        if isempty(iend)==0;
            [s,is]=sort(abs(iend-lhs(ii-1,2)));   %for the lhs
            if abs(iend(is(1))+1-iFBragg(1))<=N_max; lhs(ii,2)=iend(is(1))+1; end
            %fix this estimate if the whole FOL is on one side of the Bragg
            %and there is a second segment that wraps thru bragg that I'm missing
              if iend(is(1))-iFBragg(1) < 0  %it should be odd to have the ending value be on the wrong side of the bragg line
                if length(is)-is(1)>0 %there is an end to the right of this one
                    if iend(is(1)+1)-iFBragg(1) > 0 & abs(iend(is(1)+1)-iFBragg(1))<N_max %the iend to the right of the shortest path is on the correct side of the Bragg line
                        lhs(ii,2)=iend(is(1)+1);   %use this one instead
                    end
                end
              end
             
            [s,is]=sort(abs(iend-rhs(ii-1,2)));   %for the rhs
            if abs(iend(is(1))+1-iFBragg(2))<=N_max; rhs(ii,2)=iend(is(1))+1; end
            %fix this estimate if the whole FOL is on one side of the Bragg
            %and there is a second segment that wraps thru bragg that I'm missing
            if iend(is(1))-iFBragg(2) < 0;  %it should be odd to have the ending value be on the wrong side of the bragg line
                if length(is)-is(1)>0 %there is an end to the right of this one
                    if iend(is(1)+1)-iFBragg(2) > 0 %the iend to the right of the shortest path is on the correct side of the Bragg line
                        rhs(ii,2)=iend(is(1)+1);   %use this one instead
                    end
                end
             end

        end        
    end  
end


end

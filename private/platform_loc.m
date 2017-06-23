function hp = platform_loc(swch)
% PLATFORM_LOC.M
% hp=platform_loc(swch)
% Adds platform locations to any plot
% Called by sbchan_map.m
%
% optional input swch=0 turns of labels.
% Output's handle of platform objects.

% made, overly complexly, by Brian Emery

% input illogic:
if nargin==0
    swch=0; % default no labels
elseif nargin==1 & swch==1
    swch=1;
elseif nargin==1 & swch==0
    swch=0;    
end

hp=[];

%Site	Latitude	Longitude
site=struct('name',{'Gina','Gail','Gilda','Grace','Habitat','Hogan','Houchin','Henry',...
                    'Hillhouse','A','B','C','Holly','Hondo','Harmony','Heritage','Hermosa', ...
                    'Harvest','Hidalgo','Irene'}, ...
            'lat',{34.1175,34.125,34.18216667,34.1795,34.28666667,34.33783333,34.335,34.33333333, ...
                    34.33133333,34.332,34.33216667,34.333,34.38966667,34.39083333,34.377,34.35,34.4555, ...
                    34.46916667,34.495,34.61033333}, ...
            'lon',{-119.2763333,-119.4003333,-119.4178333,-119.4686667,-119.5883333,-119.5413333,-119.5521667,-119.5603333, ...
                    -119.6033333,-119.6125,-119.6216667,-119.6308333,-119.9055,-120.1205,-120.1675,-120.2791667,-120.6463333, ...
                    -120.6808333,-120.7021667,-120.7233333});
			
for i=1:20
    hi=plot(site(i).lon,site(i).lat,'ks');,set(hi,'Clipping','on','MarkerSize',3,'MarkerFaceColor',[.5 .5 .5]);
    hp=[hp hi];
    if swch
    text(site(i).lon,site(i).lat+.015,zeros(size(site(i).lon)),site(i).name,'fontsize',8,'Clipping','on');
    end
end

    
% original info:    
%    No.	Site	Latitude	Longitude
%1	Gina	34.1175	-119.2763333
%2	Gail	34.125	-119.4003333
%3	Gilda	34.18216667	-119.4178333
%4	Grace	34.1795	-119.4686667
%5	Habitat	34.28666667	-119.5883333
%8	Hogan	34.33783333	-119.5413333
%9	Houchin	34.335	-119.5521667
%10	Henry	34.33333333	-119.5603333
%11	Hillhouse	34.33133333	-119.6033333
%12	A	34.332	-119.6125
%13	B	34.33216667	-119.6216667
%14	C	34.333	-119.6308333
%17	Holly	34.38966667	-119.9055
%18	Hondo	34.39083333	-120.1205
%19	Harmony	34.377	-120.1675
%20	Heritage	34.35	-120.2791667
%21	Hermosa	34.4555	-120.6463333
%22	Harvest	34.46916667	-120.6808333
%23	Hidalgo	34.495	-120.7021667
%24	Irene	34.61033333	-120.7233333

end
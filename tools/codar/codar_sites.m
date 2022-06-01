function loc = codar_sites(sites)
% CODAR_SITES.M - structure containing site locations
% locn =codar_sites
%
% Outputs the long and lat of the codar sites, and the
% midpoint between ptc and cop, which the grid is based on,
% in a structure called loc.
%
% This mfile essentially makes this information
% easy to change, or add to, without modifying 
% the many mfiles that use the information
%
% now can limit the outputs if given an input cell:
% loc = codar_sites({'cop','central','ssd'})
%
% loc = 

% New Verison: Structures!

% Brian Emery 5May2008

%  Define the locations of the radar sites, and the central lat/lon
%  to convert all the radial files and total grid data to.
loc=struct();



% ----
% UCSB SITES

% This is from field notes dated 16Dec1997. Its the same as those recorded
% in "Settings at Sites" in the back of Field notebook Vol.1, but that one
% is from an unknown source.
loc.ptc = [(-120-(28.299/60)) 34+(26.897/60)]; 

% See page 1 of Field notebook Vol. 1. There is no record of where the
% position (-119.8767 34.4100) is from. It's about 300m inland from this one:
loc.cop = [-119-(52.7239/60) 34+(24.4564/60)];  

% The position below is very close to this: (-120-(4/60)-(37/3600)),(34+(27/60)+(41/3600))
% which is recorded in "Settings at Sites" in the back of Field notebook Vol.1,
% but I have no record of it's measurement.
loc.rfg= [-120.0767 34.4612];

% From Notes dated 23Nov96. This [(-120-(38.921/60)) (34+(34.608/60))] is about
% 160m from location given on pg 95 of Field notebook Vol.1, which is stated
% below
loc.arg=[(-120-(39.027/60)) (34+(34.605/60))] ;

%from 30nov99 measurement (pg 148 of Field notebook Vol.1)
loc.fbk=[(-120-(37.260/60)) (34+(52.170/60))]; 

% Oxnard's Mandalay Generating Station
loc.mgs=[-119-(15.083/60) 34+(12.304/60)];

% summerland site May-Aug 2004
loc.smr=[-119.60393 34.420433];

% summerland site march 2005, summerland sanitary district
loc.ssd=[-119-(35.765/60) 34+(25.140/60)];

% SCI NAVY SITE (see SCI_GPS.txt)
% -119.6312284, 33.9949065 % this must be old
% this from header file (March 2007 measurement)
loc.sci=[-119.63116667 33.9948167];

% POINT MUGU!
loc.ptm = [-119.10736 34.09612]; % GPS of rcv measured 27Aug2008

% Nicholas Canyon
loc.nic = [-118-(54.918/60) 34+(02.546/60)];

% Gaviota 
loc.gta=[-120-(12.46/60) 34+(28.03/60)];

% SNI
% 33�16.81583'N,119�30.6705'W
loc.sni=[-119.5224500 33.2805000];

% Santa Rosa?
loc.sri=[-120.1157000 33.8918000];

loc.trl = [-120.2322   34.4711];

% Hollister Ranch Gate
loc.hol = [-120-(14.941/60) 34+(28.134/60)];  

% Enable 4 letter site codes too ...
fn = fieldnames(loc);

for i = 1:numel(fn)
    loc.([fn{i} '1']) = loc.(fn{i});
end




% ----
% SLO SITES
loc.dcsr = [-120-(50.764/60) 35+(12.148/60)];
loc.dclr = [-120.8628 35.2172];
loc.luis = [-120.7564 35.1601];
loc.ragg = [-121.3362500 35.7873610];
loc.agl1 = [-120-(38.945/60) 34+(34.613/60)];
loc.estr = [-120.9776 35.4597];


% ----
% for tests ..
loc.clhb = [-65.6546500 43.4559833  ];
loc.fort = [-122.5012333 37.7120833 ];

% RFG Patch
loc.Rfg1=loc.rfg;



% EXTRAS

loc.sdsc=[-118.4869000 32.9175500];


% SNI LR
% approx from google earth
%  33 deg, 13' 53.39'' N
% 119 deg, 28' 42.77'' W
loc.snlr = [-119-(26/60)-(42.77/3600) 33+(13/60)+(53.39/3600)];

% ANACAPA
% est, using ginput
loc.acap = [-119.3835   34.0104];

% Enable upper case

fn = fieldnames(loc);

for i = 1:numel(fn)
    loc.(upper(fn{i})) = loc.(fn{i});
end



% Probably obsolete:
% define the mid point between the sites. Use the old positions to
% Keep the HF total grid from changing.
loc.central = [(-120.4637+(-119.8767))/2  (34.4543+34.4100)/2]; 





% --------------------------------------------------------- 
% SCCOOS AND OTHER SITES
% --------------------------------------------------------- 

% Add Other SCCOOS sites here, to be used as needed, plus comming soon and
% other future sites ...

% % GAVIOTA
% loc.gta1=[-120.2077000 34.4672000]; % from sccoos site - verify
% % Santa Rosa?
% loc.sri1=[-120.1157000 33.8918000];
% % ormand bch:
% loc.orb1=[-119-(10.233/60) 34+(07.632)/60 ];

% 25's:
% SIO sites
loc.SCDB = [-118.7336600 34.0332000];
loc.SCDH = [-118.4423500 33.9432000];
loc.SCPF = [-118.2945000 33.7048000];
loc.SCNB = [-117.9312600 33.6059700];
loc.SCTB = [-118.3914 33.8117];

loc.SDDP = [-117.7147000 33.4602000];
loc.SDDM = [-117.5956000 33.3884000];
loc.SDSE = [-117.2861000 33.0245000];
loc.SDWW = [-117.2474990 32.6791660];
loc.SDPL = [-117.2396 32.6658];
loc.SDCI = [-117.2437 32.4141];
loc.SDBP = [-117.1223 32.5359];

loc.UABC = [-117.0758000 32.3761000];

% 13's:
loc.ESTR=[-120.9775000 35.4597500];
loc.ORB1=[-119.1666700 34.1333300];
loc.SCCI=[-118.4781833 33.4468500];

% long range:
loc.DCLR=[-120.8626670 35.2174720];
loc.RAGG=[-121.3362500 35.7873610];

loc.gpt1=[-119.8445000 34.4048000];
loc.sri2=[-120.1157000 33.8918000];

loc.SDSB=[-119.0378950 33.4760280]; % SBI
loc.SDSC=[-118.4869000 32.9175500];
loc.SDSL=[-117.2525000 32.8706000];
loc.LOMA=[-117.2439530 32.6664040];


% WHOI
% retrieved from site's config radar_pattern.m ...
loc.NWTP = [ -(70+6/60+24.85/3600)  41+14/60+30.96/3600];
loc.METS = [ -70.5272500 41.3498000 ];
loc.SQUB = [ -70.7675833 41.3065667 ];
loc.NANT = [ -69.9750 41.2506 ];
loc.MVCO = [ -70.5267667 41.3497833 ];
loc.LPWR = [-(70+38.41/60) 41+20.904/60 ];


% THIS MUST BE LAST
% Now limit outputs if there is an input
if nargin ==1
    loc = rmfield(loc, setdiff(fieldnames(loc),sites) );  
    
end






return
% --------------------------------------------------------- 
% POSSIBLE NEW SITES
%  --------------------------------------------------------- 
% also a record of old sites and site surveys

% from survey Nov 2001: End of Carp Pier 34 23.053, 119 30.475
% car=[-119-(30.475/60)	34+(23.053/60)];
%
% Rincon Island (this from NOAA Chart, which quotes the location)
% plot(-119-(26.6/60),34+(20.8/60),'g*')
%rci=[-119-(26.6/60) 34+(20.8/60)];

% Santa Barbara Light, pg 107 in field notebook
% sbl=[(-119-43.368/60) (34+23.758/60)];
%
% Gaviota Processing Facility (roughly)
% gav=[-120-(12.46/60) 34+(28.03/60)];     

% for reference: Ellison Hall
% ellison=[(-119-50.697/60) (34+24.914/60)];

% ormand beach (and fence line)
% plot([-119-(10.233/60) -119-(10.184/60)],[34+(07.632)/60 34+(07.697)/60],'g*')


% PRE 11Apr00 CODE:
%  Define the locations of the radar sites, and the central lat/lon
%  to convert all the radial files and total grid data to.
ptc = [-120.4637 34.4543]; 
cop = [-119.8767 34.4100];    
rfg= [-120.0767 34.4612];

% From ? ~160m from location given on pg 95 of Field notebook Vol.1
arg= [(-120-(38.921/60)) (34+(34.608/60))];

%from 30nov99 measurement (pg 148 of Field notebook Vol.1)
fbk=[(-120-(37.260/60)) (34+(52.170/60))]; 

% define the mid point between the sites
central = [(ptc(1)+cop(1))/2  (ptc(2)+cop(2))/2];  
end

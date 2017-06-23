function [Vr,Bear,Range] = compute_range_cell_velocity(bear,nRC,TUV,siteloc)
% COMPUTE RANGE CELL VELOCITY - output the radial velocity vs range cell
% 
% -- this is sort of obsolete !!! ----
%
% Need to think about how Codar is doign this vs how I want to do it. The
% signal should be ~ the sum of many small emitters? (see draft nsf
% proposal)
%

% TO DO
% need to get CFG from the simlation so i know how many and what size the 
% radial bins are. Then finish compute_rc_vel_from_tuv
% 
% What is the output of this function? probably need bearing info with it too?
%


if nargin < 3 || isempty(TUV)
    
    % TUV not given or empty, use codar's method
    % here, nRC is the range cell number 
    Vr = compute_rc_vel(bear,nRC); [Bear,Range] = deal([]);
    
else
   
    % Do all range cells? Need a site location
    [Vr,Bear,Range] = tuv_to_radcomp(TUV,siteloc,nRC);
    
    
end



end

function Vr = compute_rc_vel(bear,nRC)
% CODARS METHOD (from ComputeRCPat.m)

Nth = length(bear);

%// Compute the radial current patterns
%
%
% global simT1r1;			%//1st RC, Current fundamental direction
% global simT1r2;			%//1st RC, Current 2nd harmonic direction
%
% global simV2m0;			%//Last RC, ...
% global simV2m1;			%//Last RC, ...
% global simV2m2;			%//Last RC, ...
%
% global simT2r1;			%//Last RC, ...
% global simT2r2;			%//Last RC, ...		
% 
% global simVr1;			%//Radial current patterns for first range cell
% global simVr2;			%//Radial current patterns for last range cell



%//Default current parms
% global simV1m0;			%//1st RC, Current constant component
% global simV1m1;			%//1st RC, Current 1st harmonic strength
% global simV1m2;			%//1st RC, Current 2nd harmonic strength
simV1m0 = 0.0;
simV1m1 = 40.0;
simV1m2 = 0.0;

% simT1r1 = 0.0;
% simT1r2 = 0.0;

simV2m0 = 0.0;
simV2m1 = 40.0;
simV2m2 = 0.0;

simT2r1 = 0.0;
simT2r2 = 0.0;

%//Current standard deviation about mean direction (cm/s)
simCStd = 1.0; 


%// Use parameter equation to compute radial velocities
% if simOpt(7) == 0
%// The mean current pattern for the first range cell

% CURRENT PATTERN
% Create a flow from the east, moving west, parallel to the coast

% This affects how the currents get computed
simTh = (bear - 90)*pi/180;

% tr1 = simT1r1*pi/180;
tr2 = 0 ; %simT1r2*pi/180;
simVr1 = simV1m0 + simV1m1*cos(simTh) + simV1m2*cos(2*(simTh-tr2));

%// The mean current pattern for the last range cell
tr1 = simT2r1*pi/180;
tr2 = simT2r2*pi/180;
simVr2 = simV2m0 + simV2m1*cos(simTh) + simV2m2*cos(2*(simTh-tr2));


% % THIS MOVED TO COMPUTE_SPECTRA
% %// Add current noise
% cnoise = randn(Nth,1)*simCStd;
% 
% isNoise(2) = 1; % < --- move this to CFG
% 
% if isNoise(2) == 1
% 	%//Randomize the current (interpolated for range cell).  <-------- RANGE CELL INTERP
% 	Vr = (((31-nRC)/30)*simVr1 + ((nRC-1)/30)*simVr2) + cnoise;
% else
% 	%//Interpolate the current for desired range cell
	Vr = ((31-nRC)/30)*simVr1 + ((nRC-1)/30)*simVr2;
% end



% %// Use user supplied function to compute radial velocities
% else
%     
%     keyboard
%     
% 	%// Convert 'coast-normal' bearings to 'North' bearings (degrees)
% 	bearings = round(simTh .* (180/pi) + simCoastN);
% 
% 	%// Compute the first range cell of radial velocities
% 	dist = simZRC + (simFirstRC - 1) * simDRC;
% 	[simVr1,err] = CSSimRadVel( dist, bearings );
% 
% 	%// Compute the last range cell of radial velocities
% 	dist = simZRC + (simLastRC - 1) * simDRC;
% 	[simVr2,err] = CSSimRadVel( dist, bearings );
% end

% % EXPERIMENTS
% Vr = Vr + 40;
% 
% Vr = 50:-0.5:-50;
% Vr = Vr(1:length(bear))';
% 

% UNIFORM FLOW-ISH
% Vr = 30*sind(bear); 
% Vr = (30+randn(size(bear)).*3).*sind(bear);



end


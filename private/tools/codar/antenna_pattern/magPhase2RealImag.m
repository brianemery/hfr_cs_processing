function APM = magPhase2RealImag(APM)
% MAG PHASE 2 REAL IMAG - APM real and imaginary from mag and phase
%
% APM = magPhase2RealImag(APM)
% APM must contain the fields 'A13M','A23M','A13P','A23P',
% where phases are in degrees.
%
% Optionally does s 'A12M','A12P',
%
% See also realImag2MagPhase, abs, angle.

% Copyright (C) 2010 Brian M. Emery
% verified method with mag_phase_calc_experiments.m

if strcmp('--t',APM), test_case, return, end

field_check(APM,{'A13M','A23M','A13P','A23P'})

% COMPUTE REAL AND IMAGINARY COMPONENTS
% r cos(p) + r i sin(p)
% real = r cos(p)
% imag = r sin(p)
APM.A13R = APM.A13M .* cosd(APM.A13P);
APM.A13I = APM.A13M .* sind(APM.A13P);

APM.A23R = APM.A23M .* cosd(APM.A23P);
APM.A23I = APM.A23M .* sind(APM.A23P);

% Optional cross terms
if all(isfield(APM,{'A12M','A12P'}))
    
    APM.A12R = APM.A12M .* cosd(APM.A12P);
    APM.A12I = APM.A12M .* sind(APM.A12P);
    
end


end

function test_case
% TEST CASE
% originally verified method with mag_phase_calc_experiments.m
% make sure this works by testing on a 360 deg ideal pattern

% load ideal 360 deg pattern
APM = make_ideal_pattern(0);

% plot
plot_apm_polar(APM)

% recompute the real and imag components
APM = magPhase2RealImag(APM);

% re-plot
figure
plot_apm_polar(APM)


keyboard


end
function APM = realImag2MagPhase(APM)
% REAL IMAG 2 MAG PHASE - APM mag phase from real and imag components
%
% APM = realImag2MagPhase(APM)
% APM must contain the fields 'A13R','A23R','A13I','A23I'. Outputs phases
% in degrees (Magnitude in ? Volts I think )
%
% See also magPhase2RealImag, abs, angle, volts2dbm

% Copyright (C) 2010 Brian M. Emery
% verified method with mag_phase_calc_experiments.m


field_check(APM,{'A13R','A23R','A13I','A23I'})

% COMPUTE Magnitudes
APM.A13M = sqrt( (APM.A13R.^2) + (APM.A13I.^2) );
APM.A23M = sqrt( (APM.A23R.^2) + (APM.A23I.^2) );

% COMPUTE Phases
APM.A13P = atan2(APM.A13I,APM.A13R).*180/pi;
APM.A23P = atan2(APM.A23I,APM.A23R).*180/pi;

end
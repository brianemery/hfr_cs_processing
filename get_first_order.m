function [tf,SNR] = get_first_order(CS,minsnr)
% GET FIRST ORDER - new unfinished potentially simpler version
% [tf,SNR] = get_first_order(CS,minsnr)
%
% 2014 IDEAS
% Use logical indexing, turn this into one big logical statement, output a
% logical firstOrder matrix, same size as CS data, 1 if first order. Seems 
% like this could be greatly simplified - but then it always does! Convert
% everything to same size matrix: velocity, smoothed CS, factor stuff ...
%
% "Could probably re-implement this as a logical/find type operation"
% TO DO
% - just be one big logical statement
% - do everything in dB from the outset
% - apply to all range cells (test for CS formatting)
% - Then: screen by max current, snr, 
% - use my SNR code
%
% see subfunction to bearing_uncertainties.m


% % DEFAULT SETTINGS
% header.flim        = 15.0; %15.85;   % these are in Volts^2 (they get converted to dBm) (Default 15.0)
% header.currmax     = 100;    % cm/s
% header.fdown       = 7.5; %5.01; %    % these are in Volts^2
% header.noiseFactor = 4.0; %10.0; %    % " " " " in Volts^2
% % number of points to smooth for finding nulls (looks like what is actually
% % used is this number *2+1! (see SpectraPlotterMap)
% header.nsm = 4; 

if nargin < 2
    minsnr = 7;
end


% dbm everything
CS = cs_volts2dbm(CS);

% Get SNR struct
SNR = get_SNR(CS);

if ~isfield(CS,'Vrad')
    [~,CS.Vrad,~] = getVelocities(CS.Header);
end

% Get matrix of (ocean current) velocities
% Note rows are fbins, and columns are range cells typcially
Vrad = CS.Vrad(:)*ones(1,size(CS.antenna3Self,2));


% Create the logical array
tf = (abs(Vrad) <= 100 & ...
     SNR.antenna3Self > minsnr & (SNR.antenna13CrossSp > minsnr | SNR.antenna23CrossSp > minsnr));



% % Check plot
% figure
% plot(CS.freqs,CS.antenna3Self(:,1),'-r.')
% hold on
% plot(CS.freqs(tf),CS.antenna3Self(tf,1),'*b')
% 
% keyboard
% 


end

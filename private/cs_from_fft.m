function CS = cs_from_fft(V1,V2,V3,CFG)
% FORM THE CROSS SPECTRA - from FFTs
%
% INPUTS
% V1,V2,V3 - Doppler spectra, ie the output of the FFT performed on the
%            Range data 
% CFG      - Optional configuration specifying some of the outputs


% NOTES
% "In most practical situations, this level of information is not 
% available, however, and one must resort to estimating a given signal's
% power spectral density. One very straightforward approach is to take 
% the magnitude squared of its Fourier transform (or, perhaps, take the 
% squared magnitude of several short-time Fourier transforms and average
% them) as the estimate of the PSD."
%
% So the Auto spectra are the PSD!

% set defaults
if nargin < 4
    CFG.output_ts = false;
    CFG.output_fft = false;
end

% use subfunction to make a header
Header = cs_header_struct;


if nargin == 4
    
    % Reconcile range and CS header data
    Header.freqMHz     = CFG.txfreq;
    Header.SwBWkHz     = CFG.bw;
    Header.SwRfreqHz   = CFG.srf;
    Header.nRangeCells = CFG.nRC;
    Header.fftLength   = CFG.nfft;
    
end





% % Decorrelation experiment code
% %
% % Also sum over all the velocities and divide by a factor to get the
% % overall signal level right. This should cancel out the summing
% n = size(V1,1);
% 
% 
% %// Compute the first set of cross spectra 
% antenna1Self = sum( abs(V1).^2 ).'./n;
% antenna2Self = sum( abs(V2).^2 ).'./n;
% antenna3Self = sum( abs(V3).^2 ).'./n;
% antenna12CrossSp = sum( V1.*conj(V2) ).'./n;
% antenna13CrossSp = sum( V1.*conj(V3) ).'./n;
% antenna23CrossSp = sum( V2.*conj(V3) ).'./n;

antenna1Self = abs(V1).^2 ;
antenna2Self = abs(V2).^2 ;
antenna3Self = abs(V3).^2 ;
antenna12CrossSp = V1.*conj(V2) ;
antenna13CrossSp = V1.*conj(V3) ;
antenna23CrossSp = V2.*conj(V3)  ;


% Quality array from zero to one in value.
[ spectraQualNum ] =  deal(zeros(512,1)); 


 
 



% 
% antenna1Self = S11;
% antenna2Self = S22;
% antenna3Self = S33;
% antenna12CrossSp = S12;
% antenna13CrossSp = S13;
% antenna23CrossSp = S23;

rdCSErr = 0;
FileName = '';

% Pack data into structure
vars ={'Header','antenna1Self','antenna2Self','antenna3Self', ...
       'antenna12CrossSp','antenna13CrossSp','antenna23CrossSp', ...
       'spectraQualNum','rdCSErr','FileName','CFG'};

if CFG.output_ts
     
    vars{end+1} = 'v1'; 
    vars{end+1} = 'v2';
    vars{end+1} = 'v3';
    
end

if CFG.output_fft
       
    vars{end+1} = 'V1';  
    vars{end+1} = 'V2';
    vars{end+1} = 'V3';
    
end
      
CS = struct_pack(vars);

% I don't want to do this just yet
% CS = cs_volts2dbm(CS);

[CS.freqs,CS.Vrad,fb] = getVelocities(CS.Header);


 %[nBar,ix] = get_noise_level(CS.antenna3Self,CS.freqs,fb);
 
 
 

end



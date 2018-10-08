function CS = cs_from_fft(fftn,CFG) %V1,V2,V3,CFG)
% FORM THE CROSS SPECTRA - from FFTs
%
% INPUTS
% fftn     - The N rows by length(FFT) columns 
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


% check for test case
if strcmp('--t',fftn), test_case, return, end


% set defaults
if nargin < 2
    CFG.output_ts = false;
    CFG.output_fft = false;
end

% get number of antennas
m = size(fftn,1);

if m == 3
    
    % Get struct setup
    CS = cs_struct;
    
else
    CS = cs_struct(1,m);
    
end


% use subfunction to make a header
CS.Header = cs_header_struct;


if nargin == 2
    
    % Reconcile range and CS header data
    CS.Header.freqMHz     = CFG.txfreq;
    CS.Header.SwBWkHz     = CFG.bw;
    CS.Header.SwRfreqHz   = CFG.srf;
    CS.Header.nRangeCells = CFG.nRC;
    CS.Header.fftLength   = CFG.nfft;
    
end




if m == 3
    
    % use previous method and formats ... wants columns
    CS = old_way(fftn(1,:).',fftn(2,:).',fftn(3,:).',CS);
    
else
    
    % new way - arbitrary arrays
    CS = new_way(fftn,CS);
    
    
end


% add couple things

% Quality array from zero to one in value.
[ CS.spectraQualNum ] =  deal(zeros(size(fftn,2),1)); 

CS.rdCSErr = 0;
CS.FileName = '';

[CS.freqs,CS.Vrad,~] = getVelocities(CS.Header);


if CFG.output_fft
       
   % need to fix this .. or punt!
   %  CS.V1 = V1;
    CS.fftn = fftn;
    
end




end

function CS = old_way(V1,V2,V3,CS)

CS.antenna1Self = abs(V1).^2 ;
CS.antenna2Self = abs(V2).^2 ;
CS.antenna3Self = abs(V3).^2 ;
CS.antenna12CrossSp = V1.*conj(V2) ;
CS.antenna13CrossSp = V1.*conj(V3) ;
CS.antenna23CrossSp = V2.*conj(V3)  ;


end

function CS = new_way(fftn,CS)
% make the field names? See cs_struct.m ... maybe need a tool for this

% get number of antennas
m = size(fftn,1);

% Apply transpose to row shape (CS are nfft x nRC) typically
fftn = fftn.'; 


% get field names (could also read these?), and corresponding row,col
% indecies
[fn,I,J] = cs_make_field_names(m);


for i = 1:length(I)
        
        % ... ok, very slight differences from old way unless I do this
        if I(i) == J(i)
            
            CS.(fn{i}) = abs(fftn(:,I(i))).^2 ;
            
        else
            
            % now calc auto and cross products
            CS.(fn{i}) = fftn(:,I(i)) .* conj(fftn(:,J(i))) ;
            
        end
end


end


function test_case
% TEST CASE
%
% Test new method vs old method functionality
%
% data for this created with radar_simulation.m test
%



% load data
load  /m_files/test_data/cs_from_fft.mat


% Create data using old way
CSo = old_way(fftn(1,:).',fftn(2,:).',fftn(3,:).',CS);

fno = cs_fieldnames(CSo);

% Create data using new way
CSn = new_way(fftn,CS);


% compare
% fno = 
% 
%     'antenna1Self'
%     'antenna2Self'
%     'antenna3Self'
%     'antenna12CrossSp'
%     'antenna13CrossSp'
%     'antenna23CrossSp'
% 

fnn = {'a11','a22','a33','a12','a13','a23'};

disp('cs_from_fft.m ...')

for i = 1:length(fnn)
    compare_results(CSn.(fnn{i}),CSo.(fno{i}))
end

keyboard

end


function compare_results(old,new)
    
    
    if isequal(old,new)
        disp(' ... ok')
    else
        disp(' ... NOT OK')
        keyboard
    end



end
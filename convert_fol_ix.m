function peakIdx = convert_fol_ix(FOregi)
% CONVERT FOL IX - convert output of ImageFOLs to cell based format
%
% INTPUT
% FOregi -- The range/Doppler velocity indices of all points marked as
%           First Order (FO) in a cross spectra data matrix (sized 
%           number of range cells by FFT length). 
%
% OUTPUT
% peakIdx - A cell array (with numel = (# range cells)), each containing the
%           row index points within the First Order peaks.
%
% OUTPUT
% peakIdx - one range cell per cell index, each containing an array of the
% indecies of the first order 
%

% Copyright (C) 2019 Brian Emery 

rng = unique(FOregi(:,1));

peakIdx = cell(length(rng),1);

for i = 1:length(rng)
   
    peakIdx{i} = FOregi(FOregi==rng(i),2)';
    
end


end
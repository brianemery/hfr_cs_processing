function H = cs_header_struct
% MAKE CS HEADER - Make the header for a  CS file
% Header = cs_header_struct
% 
% SEE WRITE CS HEADER.M
%
% Note that the file name appears to be the start time of the data in the
% file. 

%//For dateSec, 695422 == Jan 1, 1904, 4294967296 == 2^32
%//dateSec = fix((daysAD - 695422)*(24*60*60));	%// int32
H.dateTimeSec04 =  (now - 695422)*(24*60*60) - 4294967296;


H.kindOfData =  1;           % 1 is 'raw', 2 is averaged
H.siteStr =  '';
H.averagingTimeMin =  4;
H.freqMHz =  13.4900;
H.SwRfreqHz =  2;
H.SwBWkHz =  100.7080;
H.sweepUp =  0;
H.nRangeCells =  56;
H.firstRangeCell =  0;
H.distToFirstRangeCell =  1.4895;
H.dTHex =  [];
H.deleteRawSpectra =  0;
H.overrideRemHeader =  0;
H.fftLength = 512;

% 
% %// Collect header data
% version = 4;           							%// int16
% dateSec = (daysAD - 695422)*(24*60*60) - 4294967296;
% bytesLeft1 = 62;         						%// int32
% if simOpt(5) == 1
% 	dataKind = 1;								%// int16 (Raw CS)
% else
% 	dataKind = 2;								%// int16 (Averaged CS)
% end
% bytesLeft2 = 56;								%// int32
% site = simSite;									%//  char, 4 bytes
% bytesLeft3 = 48;								%// int32
% if simOpt(5) == 1
% 	averagingTime = 0;							%// int32 (raw CS)
% else
% 	averagingTime = avgTime;					%// int32
% end
% deleteAfter = 0;								%// int32
% overrideHeader = 0;								%// int32
% frequency = txfreq;								%// float32
% sweepRepFreq = simSRF;							%// float32
% sweepBandwidth = simSBW;						%// float32  %//(=150.0/deltaRC)
% sweepIsUp = 0;									%// int32
% sizeOfFFT = 512;								%// int32
% numberOfRCs = nRCs;								%// int32
% firstRC = simFirstRC;							%// int32
% distToFirst = firstRC*150.0/sweepBandwidth;		%// float32
% bytesLeft4 = 0;					
% 
% fn = '/m_files/test_data/getDopplerVelocities/CSQ_cop1_08_12_06_200428.cs';
% dat = cs_read(fn,':',1);
% 
% Header = dat.Header;
% Header.siteStr = 'SCrz';
% Header.freqMHz = 13.45; %txfreq;
% Header.SwBWkHz = 100;
% Header.nRangeCells = 1;

end

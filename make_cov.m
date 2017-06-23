function C = make_cov(CS,fbin,rdx)
% MAKE COV - make covariance matrix from CS data
%
% From doa_on_range_cell.m for example

% Copyright (C) 2016 Brian Emery

if nargin < 3
    rdx = 1;
end

% fbin = find(CS.Vrad > Vr -2 & CS.Vrad < Vr +2);
% 
% H = cs_plot(CS,1); hold on
% LS = line_style_groups;
% h = cs_plot(CS,1,LS(1),H,fbin);
% 
% keyboard

% fbin = 399; % Need to fix the off by one aspect of the velocity

% build covariance matrix - These are the averaged complex voltages
% < ViVj* > with units volts^2
C(1,1) = CS.antenna1Self(fbin,rdx);
C(1,2) = CS.antenna12CrossSp(fbin,rdx);
C(1,3) = CS.antenna13CrossSp(fbin,rdx);

C(2,1) = conj(CS.antenna12CrossSp(fbin,rdx));
C(2,2) = CS.antenna2Self(fbin,rdx);
C(2,3) = CS.antenna23CrossSp(fbin,rdx);

C(3,1) = conj(CS.antenna13CrossSp(fbin,rdx));
C(3,2) = conj(CS.antenna23CrossSp(fbin,rdx));
C(3,3) = CS.antenna3Self(fbin,rdx);






end
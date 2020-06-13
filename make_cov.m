function C = make_cov(CS,fbin,rdx)
% MAKE COV - make covariance matrix from CS data
% C = make_cov(CS,fbin,rdx)
%
% From doa_on_range_cell.m for example

% Copyright (C) 2016 Brian Emery

% check for test case
if strcmp('--t',CS), test_case, return, end


if nargin < 3
    rdx = 1;
    
end

if isfield(CS,'antenna1Self')
    % old way
    
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
    
    
else

    
    % get field names
    fo = cs_fieldnames(CS);
    
    % detect number of antennas
    % solve n(n+1) = 2*length(fn) with quadratic eqn
    m = (-1 + sqrt(1 + 8*length(fo)))/2;
    
    % make the row,col indecies
    [fn,I,J] = cs_make_field_names(m);
    
    % might as well check on this
    if ~isequal(fo,fn), disp('CS field name discrepancy'), keyboard, end
    
    
    % preallocate 
    C = NaN(m);
    
    for i = 1:numel(fn)
       
        % get row and column index from the field name
        r = I(i);
        c = J(i);
        
        % CS contains the upper triangle, lower gets conj'd        
        C(r,c) = CS.(fn{i})(fbin,rdx); 
        
        if r~=c,
            % now fill in the lower triangle and conjugate
            C(c,r) = conj(CS.(fn{i})(fbin,rdx));
        end
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

% declare! (for R2014b 'bug')
fftn = [];

% load data
load  /m_files/test_data/cs_from_fft.mat

% make the cs
CS = cs_from_fft(fftn,CFG);

% make a copy with new field names
fn = cs_fieldnames(CS);

fnn = {'a11','a22','a33','a12','a13','a23'};

CSn = cs_struct(1,3);

for i = 1:numel(fn)
    % it's essential that fnn and fn refer to right fields!
    CSn.(fnn{i}) = CS.(fn{i});
end

% choose bin 158 which is in the peak
fbin = 158;

C1 = make_cov(CS,fbin,1);

C2 = make_cov(CSn,fbin,1);

if isequal(C1,C2)
    disp(' ...passed')
else
    disp(' ... NOT passed')
    keyboard
end


end
function CS = cs_volts2dbm(CS)
% CS VOLTS SQUARED TO DBM - convert CS data to dBm
% CS = cs_volts2dbm(CS)
%
% Converts the components of a Cross Spectra file from units of volts.^2 to
% dBm using volts2dbm.m

% Copyright (C) 2011 Brian Emery
% 6 May 2011


% Test case
if strcmp(CS,'--t'), test_case, return, end


% RECURSE IF MULTI-ELEMENT
if numel(CS) > 1, 
    
    % Add meta container
    [CS(1:numel(CS)).Units] = deal('');
    
    % Add container for phases 
    [CS(1:numel(CS)).Phases] = deal(struct([]));
    
    
    for i = 1:numel(CS)
        CS(i) = cs_volts2dbm(CS(i));
    end
    
    return
end



% Detect if units have been converted to dBm
if ( ~isfield(CS,'Units') || ~strcmp('dBm',CS.Units) ) && ~isempty(CS)
    
    %  disp('Converting to dB')
    
    % CONVERT UNITS

    % Define field names we'll loop over
    fn = cs_fieldnames(CS);
    
    
    for i = 1:numel(fn)
        
        [CS.(fn{i}),CS.Phases(1).(fn{i})] = volts2dbm( CS.(fn{i}) );
                
    end
    
    
    % Add meta data
    CS.Units = 'dBm';
end


end

function test_case
% TEST CASE

% functionality test ...

fn = '/m_files/test_data/compute_apm_from_csq/CSQ_cop1_08_12_06_205124.cs';

dat = cs_read(fn);

dat(1:3) = dat;



% run test
dat = cs_volts2dbm(dat);


keyboard

if isequal(dat(1),dat(3))
    disp('cs_volts2dbm: test ... ok')
else
    disp('cs_volts2dbm: test ... NOT ok')
end


end


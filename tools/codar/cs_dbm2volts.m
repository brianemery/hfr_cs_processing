function CS = cs_dbm2volts(CS)
% CS DBM TO VOLTS SQUARED- convert CS data to volts^2
% CS = cs_dbm2volts(CS)
%
% Converts the components of a Cross Spectra file from units of dbm to
% volts.^2 to dbm using dbm2volts.m
%
% Note that the default (no Units field), typically means the data is in
% volts^2.

% Copyright (C) 2017 Brian Emery
% 26 May 2017

% TODO 
% remove phase field?


% Test case
if strcmp(CS,'--t'), test_case, return, end


% RECURSE IF MULTI-ELEMENT
if numel(CS) > 1, 
    
    for i = 1:numel(CS)
        CS(i) = cs_dbm2volts(CS(i));
    end
    
    % can do this here, but cant figure out how to do it for the
    % singleton case ...
    CS = rmfield(CS,'Phases');
    
    return
end


% Detect if units have been converted to dBm ... these have to be true if
% the units have been converted to dBm
if isfield(CS,'Units') && strcmp('dBm',CS.Units)
    
    
    % CONVERT UNITS

    % Define field names we'll loop over ( generalize a*)
    fn = cs_fieldnames(CS);
    
    
    for i = 1:numel(fn)
        
        if isfield(CS,'Phases')
            
            % use the phases if they are there
            CS.(fn{i}) = dbm2volts( CS.(fn{i}), CS.Phases.(fn{i}) );
          
            % empty it out at least
            CS.Phases = rmfield(CS.Phases,fn{i}); % disp('using phases ..')
            
        else           
            CS.(fn{i}) = dbm2volts( CS.(fn{i}) );
            
        end
    end
    
    
    % Add meta data
    CS.Units = 'volts^2';
end


end

function test_case
% TEST CASE

% functionality test ...

fn = '/m_files/test_data/compute_apm_from_csq/CSQ_cop1_08_12_06_205124.cs';

dat1 = cs_read(fn);

dat(1:3) = dat1;

dat = cs_volts2dbm(dat);


% now run test for this function, and check reversal worked

dat2 = cs_dbm2volts(dat);



% PLOT RESULT
fn = cs_fieldnames(dat2);


% this suggests it's reversible to order 1e-20
for i = 1:numel(fn)
    figure(1), hold on
    plot(real(dat1.(fn{i}))-real(dat2(1).(fn{i})) ,'.')
    figure(2), hold on
    plot(imag(dat1.(fn{i}))-imag(dat2(1).(fn{i})) ,'.')

end

keyboard

% if isequal(dat(1),dat(3))
%     disp('cs_volts2dbm: test ... ok')
% else
%     disp('cs_volts2dbm: test ... NOT ok')
% end




end


function APM = loop_to_apm(S)
% LOOP TO APM - convert ctfRead of Loop files to APM standard format
% 
% Originaly used for SOO data, now for regular SS LOOP files ...
%
% SEE ALSO
% translate_loop.m ... (for use with SOO loop data )
% load_pattern_file.m 

% Updates/Breakage 7 Dec 2021 (80 yr anniversary of Pearl Harbor)
% - converting to work with SeaSonde Loop files
% - need to test with the other kind, or just use translate_loop for those

if strcmp(S.Type,'LOOP')
    
    APM = apm_struct;
    % set the APM type
    APM.Type = 'LOOP';
    
else
    % create empty APM struct
    APM = apm_struct_ext(numel(S.TIME));
    
    % set the APM type
    APM.Type = 'SOO';
end


% TRANSLATE FIELDS

% Get TimeStamp from array
APM.TimeStamp = datenum(strcat(num2str(S.DATE),num2str(S.TIME,'%06.0f')),'yyyymmddHHMMSS');

% Char fields that directly match
fn = {'SiteName','SiteOrigin','FileName'};

for i = 1:numel(fn)
    try APM.(fn{i}) = S.(fn{i})(1,:); catch, end
end

% Numeric fields that directly match
fn = {'A13M','A13P','A23M','A23P'};
 
for i = 1:numel(fn)
    APM.(fn{i}) = S.(fn{i});
end


% Input/Output fields
in  = {'TRGB','DPRV','PKRC','PKDC'};
out = {'BEAR','dopplerRadVel','rCellIdx','fbinIdx'};

for i = 1:numel(in)
    try
    APM.(out{i}) = S.(in{i});
    catch E
       disp(['No '  in{i} ' field'])
    end
end

% target radial velocity (output cm/s) 
try APM.RadVel.xbar = S.TGRV; catch, end

% stdev of target radial velocity during CSQ (cm/s)
try APM.RadVel.stdev = S.TGSD; catch, end



% APM Fields to compute
APM = magPhase2RealImag(APM);

% SNR fields  
APM.SNR.stdCodar  = [S.A1SN S.A2SN S.A3SN];
try APM.SNR.backGrnd  = [S.SBG1 S.SBG2 S.SBG3]; catch, end
try APM.SNR.localDopp = [S.SLD1 S.SLD2 S.SLD3]; catch, end
try APM.SNR.inRange   = [S.SIR1 S.SIR2 S.SIR3]; catch, end


if isfield(S,'AR3D') && isfield(S,'AR3P')
    
    % S.AR3D, S.AR3P to A33R A33I
    APM.A33R = 10.^(S.AR3D./10).* cosd(S.AR3P);
    APM.A33I = 10.^(S.AR3D./10).* sind(S.AR3P);
    
    % add these too 
    APM.A33M = 10.^(S.AR3D./10); % convert to volts
    APM.A33P = S.AR3P; % deg
    
end

APM = add_units(APM);

end


function APM = add_units(APM)


% % Meta
APM.README.BEAR_Units = 'degCWN';
APM.README.loop1Brg_Units = 'degCWN';
APM.ProcessingSteps{end+1} = 'loop_to_apm';

% Add to README:
APM.Units.BEAR = 'degCWN';

fn = {'A13M','A23M','A33M'};
for i = 1:numel(fn), APM.Units.(fn{i}) = 'volts'; end

fn = {'A13P','A23P','A33P'};
for i = 1:numel(fn), APM.Units.(fn{i}) = 'deg'; end





end


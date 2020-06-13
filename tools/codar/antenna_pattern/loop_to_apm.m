function APM = loop_to_apm(S)
% LOOP TO APM - convert ctfRead of Loop files to APM standard format
%
% * intended for use with SOO loop data *
%
% SEE ALSO
% translate_loop.m ... need to merge these together. See also
% load_pattern_file.m 

% create empty APM struct
APM = apm_struct_ext(numel(S.TIME));

% set the APM type
APM.Type = 'SOO';


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

% S.AR3D, S.AR3P to A33R A33I
APM.A33R = 10.^(S.AR3D./10).* cosd(S.AR3P);
APM.A33I = 10.^(S.AR3D./10).* sind(S.AR3P);


end


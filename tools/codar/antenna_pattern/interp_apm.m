function Ai = interp_apm(A,BEAR)
% INTERP APM - interpolate APM data to given bearings
% Ai = interp_apm(A,BEAR)
%
% Bearings in A should equal or exceed the range of bearings in BEAR
%
% OUTPUTS
% Ai = APM struct with fields same size as BEAR
%
% TODO
% change name to apm_interp.m
%
% SEE ALSO
% struct_interp.m

% Copyright (C) 2011 Brian Emery
% 
% updates 18 Feb 2022
% - column vector output, test check code

% TO DO
% % should error or something if the APM bearings dont cover range(BEAR)?
%
% NOTE 1
% at this moment the goal of this function is to expand the APM to the size
% of the BEAR input. Thus, unique doesn't matter, just lookup the table at 
% ever given bearing. DOes this goal change? 
%
% NOTE 2
% extrap must be used on line 41 when bearings are even just a bit off

% % check just in case (sorts also)
% BEAR = unique(round(BEAR*10)/10); % SEE NOTE (1) 

% Transfer and update meta data
Ai = A;
Ai.Type = [A.Type ' (Interpolated)'];
Ai.BEAR = BEAR;


% Get list of fields to interp
fn = {'A13R', 'A13I', 'A23R', 'A23I','A33R','A33I','A'}; 


% do the interpolation
for i = 1:numel(fn)
    
    if isfield(A,fn{i}) && ~isempty(A.(fn{i}))
        Ai.(fn{i}) = interp1(A.BEAR,A.(fn{i}),Ai.BEAR,'linear','extrap'); 
        
        % make sure it's column vectors
        Ai.(fn{i}) = Ai.(fn{i})(:);
        
    end  
    
end




% Recompute Mag and Phase
try
    Ai = realImag2MagPhase(Ai);
catch
end

% clear fields not (yet) used
fn = {'A13RQ','A13IQ', 'A23RQ','A23IQ'};

for i = 1:numel(fn)
    try
        Ai = rmfield(Ai,fn{i});
    catch
    end
end


% Update processing info
Ai.ProcessingSteps{end+1} = mfilename;


end

function test_case

apm_file = ['/projects/hf_winds/data/cop1/RadialConfigs/MeasPattern.txt'];

APM = load_pattern_file(apm_file);

M = interp_apm(APM,min(APM.BEAR):0.1:max(APM.BEAR));

% fn = {'A13M','A13P','A23M','A23P'};

LS = get_linestyles(8);


fn = {'A13R', 'A13I', 'A23R', 'A23I'};

hx = plot_struct(M,'BEAR',fn,LS(4:7));
hold on
hx = plot_struct(APM,'BEAR',fn,LS(1:5),hx); hold on



end
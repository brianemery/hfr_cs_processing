function Ai = interp_apm(A,BEAR)
% INTERP APM - interpolate APM data to given bearings
% Ai = interp_apm(A,BEAR)
%
% Bearings in A should equal or exceed the range of bearings in BEAR
%
% OUTPUTS
% Ai = APM struct with fields same size as BEAR

% Copyright (C) 2011 Brian Emery

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
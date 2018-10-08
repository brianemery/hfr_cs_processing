function Ai = interp_apm(A,BEAR)
% INTERP APM - interpolate APM data to given bearings
% Ai = interp_apm(A,BEAR)
%
% Bearings in A should equal or exceed the range of bearings in BEAR

% Copyright (C) 2011 Brian Emery

% TO DO
% % should error or something if the APM bearings done cover range(BEAR)?

% check just in case (sorts also)
BEAR = unique(BEAR);

% Transfer and update meta data
Ai = A;
Ai.Type = [A.Type ' (Interpolated)'];
Ai.BEAR = BEAR;


% Get list of fields to interp
fn = {'A13R', 'A13I', 'A23R', 'A23I','A33R','A33I','A'}; 


% do the interpolation
for i = 1:numel(fn)
    
    if isfield(A,fn{i}) && ~isempty(A.(fn{i}))
        Ai.(fn{i}) = interp1(A.BEAR,A.(fn{i}),Ai.BEAR,'linear'); %,'extrap');  
    end  
    
end


% Recompute Mag and Phase
Ai = realImag2MagPhase(Ai);

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
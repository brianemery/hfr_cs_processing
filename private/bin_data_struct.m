function T = bin_data_struct(S,fld,bw,bc)
% BIN DATA STRUCT - implementation of binData for structures
% T = bin_data_struct(S,fld,bw,bc)
%
% INPUT
% A structure (APM)
% the field containing bin variable
% the bin width and optionally bin centers
%
% "% x can be a matrix as long as the number of
% columns is equal to the length of y."
%
% See Also: binData, make_analysis_plots, bin_apm_struct
%
% Needs generalization, in particular, need to be able to specify rows or
% columns to bin over ... though this is typically used to bin in time so
% maybe assume it's over columns?
%
% See also, bin_radials.m (specific implementation) and it's subfunction
% which flips dimensions

% Copyright (C) 2010 Brian M. Emery

% TO DO
% make less specific to make_analysis_plots implementation

% get input fields
fn = fieldnames(S);

% Create output structure
T = S; T.Meta = struct([]);

% exclude bin variable field
fn = setdiff(fn,fld);

% loop over fields applying binning
for j = 1:numel(fn)
    
    % make sure we only do this to the right fields
    if isnumeric(S.(fn{j})) && (size(S.(fn{j}),2) == size(S.(fld),2))
        
        % Create the meta sub struct 
        T.Meta(1).(fn{j}) = struct([]);
                
        % Bin and store
        [T.(fn{j}),T.(fld),T.Meta.(fn{j})(1).N, ...
            T.Meta.(fn{j})(1).stdev, T.Meta.(fn{j})(1).stats] = binData(S.(fn{j}),S.(fld),bw,bc);
%         
%         % GET DEGRESS OF FREEDOM
%         % N is number of data points in each bin, but confidence intervals need
%         % degrees of freedom. I'm not sure that adjacent measurements from an
%         % individual ship (ie within the same CSQ peak) would be considered 
%         % independent. SO, compute degrees of freedom from N and number of ships
%         T.Meta.(fn{j}).Nd = get_dof(T.Meta.(fn{j})(1).stats.index,S.shipID);
%         
    end
end

% Optionally add meta data
if isfield(T,'ProcessingSteps')
    T.ProcessingSteps{end+1} = mfilename;
end


return

% % re-compute magnitude and phase parts to avoid bin averaging acrros
% 0-360 for example
T = compute_mag_phase(T);

return
% % CODE TO 
% % think about how to bin average cyclical quantitites, then apply scaling,
% % then compute mag and phase?
% TT = compute_mag_phase(T);
% 
% % Check plots
% figure
% subplot(2,1,1)
% check_plot(S,T,TT,'A13M')
% 
% subplot(2,1,2)
% check_plot(S,T,TT,'A13P')
% 
% figure
% subplot(2,1,1)
% check_plot(S,T,TT,'A23M')
% 
% subplot(2,1,2)
% check_plot(S,T,TT,'A23P')
% 
% keyboard

end

% -----------------------------------------------------------------------
function check_plot(S,T,TT,fld)

h(1) = plot(S.BEAR,S.(fld),'r.');
hold on
h(2) = plot(T.BEAR,T.(fld),'bo','MarkerFaceColor','b');
h(3) =plot(T.BEAR,T.(fld)/max(T.(fld)),'-co','MarkerFaceColor','c');

h(4) =plot(TT.BEAR,TT.(fld),'-mo','MarkerFaceColor','m');
h(5) =plot(TT.BEAR,TT.(fld)/max(TT.(fld)),'-mo')%,'MarkerFaceColor','m');

ylabel(fld)
xlabel('Bearing (cwN)')

legend(h,'Raw','Binned','Normalized','Recalc','Recalc Normalized')

end

% -----------------------------------------------------------------------
function  Nd = get_dof(idx,shipID)
% GET DEGRESS OF FREEDOM
% N is number of data points in each bin, but confidence intervals need
% degrees of freedom. I'm not sure that adjacent measurements from an
% individual ship (ie within the same CSQ peak) would be considered 
% independent. SO, compute degrees of freedom from N and number of ships
%
% Works with bin_data_struct.m
%
% idx is cell input of bin indexes.

Nd(1:numel(idx))=NaN;

for i = 1:numel(idx)
    
    % find the number of individual ships going into the bin average
    Nd(i) = length(unique(shipID(idx{i})));
end

end

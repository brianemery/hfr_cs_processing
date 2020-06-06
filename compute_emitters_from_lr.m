function em = compute_emitters_from_lr(LR,cut)
% COMPUTE EMITTERS FROM LR - find number of emitters given matrix of LR
% [em,tf] = compute_emitters_from_lr(LR,cut)
%
% .. this returns the first column that meets the
%    threshold ... thus minimizing the number of emitters
%
% ... also returned is a boolean that is true if the threshold was met
%     (and thus false if NOT met) ... in theory ... not done yet?
%
% Assumes:
% % get LR data to use ... could be an input? 
% LR = -2*log(S.LR.ML);
% ... this means n when LR drops below the threshold cuttoff value is
% considered the correct value of the number of emitters. (Alternatively,
% this means the LR becomes larger if the -2*log is not applied - the
% treshold is crossed going the other direction. 
% 
% OUTPUT
% em - array that gives the number of emitters for each row of a DOA matrix
%      (structured as in ?)
% tf - ** seems to be unused **
%
% SEE ALSO
% glrt.m, apply_detection.m

if nargin < 2, cut = 200; end


% create an array that gives the number of emitters for each row
em = NaN(size(LR,1),1);
%tf = false(size(LR,1),1);


% look one row at a time
for i = 1:size(LR,1)
    
    % get a column index .. this returns the first column the meets the
    % threshold ... thus minimizing the number of emitters
    try em(i) = find( LR(i,:) < cut,1); 
    
    catch
        
       em(i) = size(LR,2); % disp(['eheck ' num2str(i)]) 
       
       %tf(i) = false; % threshold NOT met
    
    end % NAN's if N=3?!
    
    % convert to absolute?
end

end


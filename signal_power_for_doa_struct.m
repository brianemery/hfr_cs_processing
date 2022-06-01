function S = signal_power_for_doa_struct(S,APM)
% COMPUTE SIGNAL POWER FOR DOA STRUCT
% S = signal_power_for_doa_struct(S)
%
% This handles the details of applying the signal_power.m calculation to
% a DOA struct, which is the pre-detection radial data format. Expects a
% field such as this for bearing:
%  S.Bear
%  ans = 
%    struct with fields:
% 
%      ML: [1901×15 double]
%      MU: [1901×15 double]
%
% Also S.OtherMatrixVars.cov: {1901×1 cell} 
%  and S.Idx 'struct of doa method detection APM indecies'
%
% SEE ALSO calc_and_store_power, sub to doa_on_range_cell.m

% 30 Jan 2021 Brian Emery

% rename data covariance matrix, still a cell array
C = S.OtherMatrixVars.cov;

% get array manifold # elements = rows
A = APM.A.';

% get doa method field names (eg ML or MU)
fn = fieldnames(S.Bear);

% create the outputs
for i = 1:numel(fn)
    S.Pwr.(fn{i}) = NaN(size(S.Bear.(fn{i})));
end

% get number of emitters searched for previously
try
    N = size(S.Idx.(fn{1}),2);
catch
    N = S.OtherMetadata(1).CFG.Nemit;
end

% get column indexing that we'll need
col = doa_column_index(N);

% loop over cov matricies, which is also the number of matrix rows
for i = 1:numel(C)
   
    % loop over DOA methods
    for j = 1:numel(fn)
    
        % loop over the number of emitters searched for
        for k = 1:N
            
            % get the apm col indecies
            ix = S.Idx.(fn{j}){i,k};
            
            if ~isempty(ix)
                
                % RUN FOR EACH METHOD
                
                % compute the power - this outputs a vector
                Pwr = signal_power(A(:,ix),C{i},S.eigValues(i,:),1); % 
                
                %  Note that music might output less than n emiters so
                % we need to account for that possibility
                pwr = NaN(1,k);
                                
                pwr(1:length(real(Pwr))) = real(Pwr);
                
                S.Pwr.(fn{j})(i,col{k}) = pwr;
                
                
                
                % compute the power - this outputs a vector
                Pwr = signal_power(A(:,ix),C{i},S.eigValues(i,:),2); % 
                
                %  Note that music might output less than n emiters so
                % we need to account for that possibility
                pwr = NaN(1,k);
                                
                pwr(1:length(real(Pwr))) = real(Pwr);
                
                S.Pwr2.(fn{j})(i,col{k}) = pwr;
                
                
                
                % compute the power - this outputs a vector
                Pwr = signal_power(A(:,ix),C{i},S.eigValues(i,:),3); % 
                
                %  Note that music might output less than n emiters so
                % we need to account for that possibility
                pwr = NaN(1,k);
                                
                pwr(1:length(real(Pwr))) = real(Pwr);
                
                S.Pwr3.(fn{j})(i,col{k}) = pwr;
                
            end
            
                  
        end
    end
end







end

function test_case

% er, not really but it could be ...

load  /projects/hf_winds/data/signal_power_testing.mat APM S


S = signal_power_for_doa_struct(S,APM)

end
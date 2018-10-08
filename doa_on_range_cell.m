function S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,n) %,snr)
% DOA ON RANGE CELL - run DOA methods on one rance cell of CS data
% [MU,ML,WM,WF,SM] = doa_on_range_cell(CS,APM,peakIdx,rdx)
%
% Typically called by doa_on_cs.m and experiment_music_vs_mle_roms.m - that
% is both by simulation code and by cs_processing code.
%
% Make sure CS have not been converted to dBm
%
% INPUTS
% CS      - Cross Spectra Structure (can be multiple range cells)
% APM     - the APM struct, mostly not used except to make the array matrix
%           and to get the bearings
% peakIdx - array of indcies of the peak to process (non-cell in this
%            function)
% rdx     - range cell index 
% K       - snapshots for music error
% n       - Search for up to n emitters
%
% OUTPUTS
% radial structs MU and ML based on MUSIC and MLE-AP, and more
%
% SEE ALSO
% doa_on_cs.m, radial_from_cs.m, cs_processing.m
%
% NOTE
% Need to incorporate it's use into doa_on_cs.m


% TO DO
% Update for new way to handle DOA structs (see doa_struct.m)
%
% lots of terrible code to deal with the variable number of outputs that
% can be returned by music. Should have handled that better with music.m
% and then reaped the benefits here ...
%
% ... I dont like this at all, should just have bear_mu, bear_ml, etc
% ...everything else is the same ...


% test?
if strcmp('--t',CS), test_case, return, end


% get array matrix for music error (SeaSondes)
A = get_array_matrix(APM);

% % matrix of single and dual bearing solutions for both methods
% % index refers to the index of the APM
% [S(1).MU,S(1).ML,S(1).WM,S(1).WF,S(1).SM] = deal(doa_struct(length(peakIdx),n)); 
% 
% % Specify both DOA method and array type
% S.MU.Type = [APM.Type ' MUSIC']; 
% S.ML.Type = [APM.Type ' MLE-AP'];
% S.WM.Type = [APM.Type ' WMUSIC'];
% S.WF.Type = [APM.Type ' WSF'];
% S.SM.Type = [APM.Type ' SML'];

S = doa_struct(length(peakIdx),n);
S.Type = [APM.Type ' DOA Struct'];

% hmm need a way to fix this in general
S.eigValues = NaN(length(peakIdx),size(APM.A,2));

mus_err = cell(n,1);

% % Define this
% fn = {'MU','ML','WM','WF','SM'};



% loop over peak indicies
for f = 1:length(peakIdx) 
    
    % Get fbin for this iteration
    fbin = peakIdx(f);
    
    % build covariance matrix     
    C = make_cov(CS,fbin,rdx);
 
    
        
    % NEED NUMBER OF SIGNALS !! ... maybe use error to guide choice
    %
    % % ... maybe do WSF here, output number of signals
    % n = mdl(C,K);
    %
    % n = aic(C,K);

    
    
    % RUN DF METHODS
    
    [doa,ftime,idx,D,V] = run_doa_calculations(A,C,APM.BEAR,n); 
    

    % Get MUSIC Error Ests
    for i = 1:n      
        mus_err{i} = music_error(A,C,APM.BEAR,K,i,idx.mu{i});
        
        % pad with NaN
        mus_err{i} = store_err(mus_err{i},i);
        
    end
     
    
    % STORE OUTPUTS
    S = update_storage(S,doa,ftime,f);
    
%     
%     
%     for i = 1:numel(fn)
%         
%         S = update_storage(S,fn{i},doa.(lower(fn{i})), ...
%                                  ftime.(lower(fn{i})),f);
%         
%         
%         
% %         S.(fn{i}) = update_storage(S.(fn{i}), ...
% %                  doa.(lower(fn{i})),ftime.(lower(fn{i})),CS.Vrad(fbin),f);
% 
%         
%     end
    
    % Save True Velocity
    S.RadVel(f,1:size(S.RadVel,2)) = CS.Vrad(fbin);
    
    % Save MUSIC things 
    S.Err(f,:) = [mus_err{:}];   
    S.eigValues(f,:) = D(:).'; 
    S.OtherMatrixVars.eigVectors{f,1} = V;
    
    % % .. could use a clean up
    % [P,tr] = run_param_test(D,V,A,mx{2});
    
    
    % Clean up for next iteration
    [mus_err] = deal(cell(n,1));
    
   
end


end

function [doa,ftime,idx,D,V] = run_doa_calculations(A,C,th,n)
% RUN DOA CALCULATIONS - from run_radar_simulation_basic.m
%
% INPUTS
% A - full array matrix
% C - data cov matrix
% n - number of signals to look for
%
% OUTPUT
% doa - structure of doa solutions, sorted and all that
% ftime - structure of doa method run times for this function call


% Expand for all possible emitters
% mx = (m*(m+1))/2;


% INIT OUTPUTS
fn = {'ml','mu','wf','sm','wm'};
for i = 1:numel(fn)
      idx(1).(fn{i}) = {};
    ftime(1).(fn{i}) = [];
    doa(1).(fn{i}) = {}; %NaN(1, (n*(n+1))/2 );
end





% COMPUTE DOAs, TRACK TIME

% MLE = CML = DML
tic
for i = 1:n
    idx(1).ml{i,1} = mle_ap(A,C,i);
end
ftime(1).ml = toc;


% MUSIC
tic
[~,idx.mu,D,V] = music(A,C,th,n);
ftime.mu = toc;

% WSF aka MODE
tic
for i = 1:n
    [idx.wf{i,1},~,~] = wsf(A,C,i);
end
ftime.wf = toc;

% SML = AML = UML
tic
for i = 1:n
    idx.sm{i,1} = sml_ap(A,C,i);
end
ftime.sm = toc;


% W-MUSIC
tic
[~,idx.wm,~,~] = music_weighted(A,C,n);
ftime.wm = toc;



% OUTPUT NAN PADDED RESULTS
% loop over fn and cells

for j = 1:numel(fn)
    for i = 1:n

        doa.(fn{j}){i,1} = store_doa(idx.(fn{j}){i},th,i);
        
    end
end


% reshape outputs
for i = 1:numel(fn)   
    doa.(fn{i}) =  [doa.(fn{i}){:}];    
end


end

function doa_mus = store_doa(idx,th,n)
% STORE DOA
% allows for possible return of 1 brg in 2 brg situation

% create storage
doa_mus = NaN(1,n);

% insert data
doa_mus(1:length(idx)) = sort(th(idx));


    
end

function err = store_err(idx,n)
% STORE DOA
% allows for possible return of 1 brg in 2 brg situation

% create storage
err = NaN(1,n);

% insert data
err(1:length(idx)) = sort(idx);


    
end

function S = update_storage(S,doa,ftime,f)
% UPDATE STORAGE
% simplify code for putting results in structs
%
%         S = update_storage(S,fn{i},doa.(lower(fn{i})), ...
%                                  ftime.(lower(fn{i})),f);

fn = fieldnames(doa);

for i = 1:numel(fn)
    S.Bear.(upper(fn{i}))(f,:) = doa.(fn{i});
    
    S.RunTime.(upper(fn{i})) = S.RunTime.(upper(fn{i})) + ftime.(fn{i});
    
end



end

function test_case
% DEV TEST

load /m_files/test_data/doa_on_range_cell.mat

S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,n);

keyboard

end


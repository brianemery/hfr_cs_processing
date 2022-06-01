function S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,CFG) %,n,dmth) %,snr)
% DOA ON RANGE CELL - run DOA methods on one range cell of CS data
% S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,n,dmth)
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
% 
% ... These in CFG struct:
% n          - Search for up to n emitters
% dmth       - (Optional) Cell struct of doa methods {'ml','mu','wf','sm','wm'}
%               defaults to {'ml','mu'}
% mus_param  - [10 5 8]; % [20 10 3] or [40 20 2];
%
% OUTPUTS
% Struct with radial sub-structs MU and ML based on MUSIC and MLE-AP, and more
% ... note that the bearing convention output is the same as that that is
% input (eg in APM) -- if it's cwN in, then it's cwN out. 
%
% SEE ALSO
% doa_on_cs.m, radial_from_cs.m, cs_processing.m, run_cs_processing*

% CHANGE LOG
% 14 Mar 2022 - removing the sort from some of the music related code to
%               keep the metrics in line with the doa data. 
%             - Added consistent use of one method to insert data into 
%               the DOA struct, including using create_doa_column_index.m  


% TO DO
% - I wonder if consolidating struct storage would be worth doing ...
%   This function could be greatly simplified by the 
%
% - make music_param an input, setable in run_cs_proc eg
%
% - Needs integration with the music-highest method ... see
%   re_run_detection_and_compare_to_drifters.m


% test?
if strcmp('--t',CS), test_case, return, end

% Enable ability to specify DOA method to use
if nargin < 6
    dmth = {'ml','mu'}; % default to just two
    % set CODAR's MUSIC Paramters
    mus_param = [10 5 8]; % [20 10 3] or [40 20 2];
else
    dmth = CFG.dmth;
    n    = CFG.Nemit;
    mus_param = CFG.mus_param;
    
end



% get array matrix for music error (SeaSondes)
A = get_array_matrix(APM); APM.A = A.';

% get data storage structure
S = doa_struct(max([length(peakIdx) 1]),n,size(A,1),dmth);
S.Type = [APM.Type ' DOA Struct'];

% document processing ... actually keep these in the file but not every
% radial struct 
% S.OtherMetadata(1).CFG = CFG;

mus_err = cell(n,1);

% compute doa column index once for the outer loop
cdx = doa_column_index(n);

% MAIN LOOP 

% loop over peak indicies
for f = 1:length(peakIdx) 
    
    % Get fbin for this iteration
    fbin = peakIdx(f);
    
    % build covariance matrix     
    C = make_cov(CS,fbin,rdx);
 
        
    % RUN DF METHODS
    % ... one est for each DOA solution ...
    [doa,ftime,idx,D,V,wdeg,doapk] = run_doa_calculations(A,C,APM.BEAR,n,dmth); 
    
     
    
    % STORE OUTPUTS
    % 'ml' and 'mu' get converted to ML and MU here
    S = update_storage(S,doa,ftime,f,idx);    
    
    % Save True Velocity
    S.RadVel(f,1:size(S.RadVel,2)) = CS.Vrad(fbin);
        
    % save this for investigations ...
    S.OtherMatrixVars.cov{f,1} = C;
    

    
    % STORE OUTPUTS FOR SPECIFIC DOA METHODS
    
    if ismember('mu',dmth)
        
        % Get MUSIC Error Ests
        for i = 1:n
            mus_err{i} = music_error(A,C,APM.BEAR,K,i,idx.mu{i});
            
            % this replced with nan_pad 
            % % pad with NaN if necessary
            % mus_err{i} = store_err(mus_err{i},i);
            
        end
        
        mus_err = nan_pad(mus_err, cdx );
        
         
        % Save MUSIC data - these should all be the right size now
        S.Err(f,:)    = mus_err;
        S.PkPwr(f,:)  = doapk;
        S.PkWdth(f,:) = wdeg;
        
        S.eigValues(f,:) = D(:).';
        S.OtherMatrixVars.eigVectors{f,1} = V;
        
        % seasonde music params - run on the dual bearing solution
        if size(C,1) == 3
            [S.Params(f,:),S.Dual(f)] = run_param_test(D,V,A,idx.mu{2},mus_param);
        end
        
    end
    
    
    if ismember('ml',dmth)
            
        % compute likelihood ratios --for all DOA methods--
        S = calc_and_store_glrt(S,idx,A,C,K,n,f);
        
    end
    
    
    % ADD SNR 
    % if the data is there
    if isfield(CS,'SNR') && isfield(CS.SNR,'antenna3Self')
        
        S.SNR(f,:) = [CS.SNR.antenna3Self(fbin,rdx) CS.SNR.antenna13CrossSp(fbin,rdx) ...
                                 CS.SNR.antenna23CrossSp(fbin,rdx)];
            
    elseif isfield(CS,'SNR') && isfield(CS.SNR,'a13')
        
        S.SNR(f,:) = [CS.SNR.a13(fbin,rdx) CS.SNR.a23(fbin,rdx) ...
                                           CS.SNR.a33(fbin,rdx)]; 
       
    end

    
    % Clean up for next iteration
    [mus_err] = deal(cell(n,1));
    
   
end


% range cell index
S.RngIdx = rdx*ones(max([length(peakIdx) 1]),1);

% update the output bearing units 
[S.RangeBearHeadUnits{2:3}] = deal(APM.Units.BEAR);

% compute signal power at the end
% ... testing use of external function, test suggests same result
% S = calc_and_store_power(S,A,n); 
S = signal_power_for_doa_struct(S,APM);


% specify which Bragg peak things are from ... true if Approaching waves,
% which are on the right, in positive Doppler frequencies
% ... friday evening cf ... sometimes empty peakIdx causes issues
if ~isempty(peakIdx)
    apprch = (peakIdx > length(CS.freqs)/2);
    S.Apprch = repmat(apprch(:),1,size(S.RadVel,2));
else
    S.Apprch = false(size(S.RadVel));
end

return

% Test code for 'Apprch' boolean
if any(~S.Apprch), keyboard, end

plot(CS.freqs,10*log10(CS.antenna3Self(:,rdx)))
hold on
plot(CS.freqs(peakIdx),10*log10(CS.antenna3Self(peakIdx,rdx)),'o')
plot(CS.freqs(peakIdx(~S.Apprch(:,1))),10*log10(CS.antenna3Self(peakIdx(~S.Apprch(:,1)),rdx)),'b*')
plot(CS.freqs(peakIdx(S.Apprch(:,1))),10*log10(CS.antenna3Self(peakIdx(S.Apprch(:,1)),rdx)),'m*')

end

function [doa,ftime,idx,D,V,wdeg,doapk] = run_doa_calculations(A,C,th,n,fn)
% RUN DOA CALCULATIONS - from run_radar_simulation_basic.m
%
% INPUTS
% A - full array matrix
% C - data cov matrix
% n - number of signals to look for
%
% OUTPUT
% doa   - structure of doa solutions, not sorted, but as a row vectors
% ftime - structure of doa method run times for this function call
% idx   - cell, not in row format since later use requires this
% D,V   - MUSIC eigenvalues and vectors
% wdeg  - vector of MUSIC peak widths as a row vector (degrees)
% doapk - vector of MUSIC DOA peak powers (db)  as a row vector

% UPDATES 
% 14 Mar 2022
% outputting from MUSIC, DOA peak power in dB and DOA half power width deg 

% Expand for all possible emitters
% mx = (m*(m+1))/2;


if nargin < 5
    % fn = {'ml','mu','wf','sm','wm'};
    fn = {'ml','mu'}; 
end


% INIT OUTPUTS
for i = 1:numel(fn)
      idx(1).(fn{i}) = {};
    ftime(1).(fn{i}) = [];
    doa(1).(fn{i}) = {}; %NaN(1, (n*(n+1))/2 );
    [D,V,wdeg,doapk] = deal([]);
end


% COMPUTE DOAs, TRACK TIME

% MUSIC
if ismember('mu',fn)
    tic
    [doafxn,idx.mu,D,V,wdeg] = music(A,C,th,n,0.05);
    for i = n:-1:1, doapk{i} = 10*log10(real(doafxn(idx.mu{i},i))); end
    ftime.mu = toc;
end


% MLE = CML = DML
if ismember('ml',fn)
    tic
    for i = 1:n
        idx(1).ml{i,1} = mle_ap(A,C,i);
    end
    ftime(1).ml = toc;
end


% WSF aka MODE
if ismember('wf',fn)
    tic
    for i = 1:n
        [idx.wf{i,1},~,~] = wsf(A,C,i);
    end
    ftime.wf = toc;
end


% SML = AML = UML
if ismember('sm',fn)
    tic
    for i = 1:n
        idx.sm{i,1} = sml_ap(A,C,i);
    end
    ftime.sm = toc;
end


% W-MUSIC
if ismember('wm',fn)
    tic
    [~,idx.wm,~,~] = music_weighted(A,C,n);
    ftime.wm = toc;
end




% OUTPUT NAN PADDED RESULTS
% loop over fn and cells

% this here to call it once
col = doa_column_index(n);

for j = 1:numel(fn)

    % get the directions from the index,  for each emitter
    for i = 1:n
        
        doa.(fn{j}){i,1} = th( idx.(fn{j}){i} );
        
        % this replaced by nan_pad (and can be removed in future)
        % doa.(fn{j}){i,1} = nan_pad_doa(,th,i);
        
    end
    
    % now nan pad
    doa.(fn{j}) = nan_pad(doa.(fn{j}),col);

end

% nan pad wdeg and doapk for music also, also output row vectors of the 
% appropriate number of columns (idx too)
if ismember('mu',fn)
    wdeg  = nan_pad(wdeg,col);
    doapk = nan_pad(doapk,col);
end

% 
% % this replaced by nan_pad (below) 
% % 
% % % reshape outputs as 1 row and the appropriate number of columns
% for i = 1:numel(fn)   
%      doa.(fn{i}) =  [doa.(fn{i}){:}];    
% end


end

function vec = nan_pad(cellin,col)

n = numel(cellin); % number of emitters searched for

vec = NaN(1,(n*(n+1))/2 ); 

for i = 1:n
    
    if ~isempty(cellin{i})

        % sometimes, the n = 2 data only has 1 thing for example
        % find those cases with this:
        % if length(cellin{i}) < i, keyboard, end
        ix = col{i};       
        vec( ix(1:length(cellin{i})) )  = cellin{i};
        
     end
end

end

function S = update_storage(S,doa,ftime,f,idx)
% UPDATE STORAGE
% simplify code for putting results in structs
%
%         S = update_storage(S,fn{i},doa.(lower(fn{i})), ...
%                                  ftime.(lower(fn{i})),f);
%
% ... this works because the doa substruct have been padded with NaN

fn = fieldnames(doa);  

for i = 1:numel(fn)
    S.Bear.(upper(fn{i}))(f,:) = doa.(fn{i});
    
    S.RunTime.(upper(fn{i})) = S.RunTime.(upper(fn{i})) + ftime.(fn{i});
    
    % Note (14 Mar 2022)
    % this indexes the cells and is correct. I want Idx to be cell to 
    % prevent NaN indecies
    S.Idx.(upper(fn{i}))(f,:) = idx.(fn{i}); 
    
end



end

function S = calc_and_store_glrt(S,idx,A,C,K,n,f)
% CALCULATE AND STORE GLRT - deal with the DOA substructs
% S = calc_and_store_glrt(S,idx,A,C,K,n);
%
% this handles the substructs and getting the output in the right format

% to allow for disabled DOA methods, get field names from the idx struct
fn = upper(fieldnames(idx)); %S.LR);


for i = 1:numel(fn)
    
    % loop over number of emitters
    for j = 1:n

%         [S.LR.(fn{i})(f,j),  S.GM.(fn{i})(f,j),  ...
%          S.Rm.(fn{i}){f,j},  S.LR2.(fn{i})(f,j),  ] 
        [S.LR.(fn{i})(f,j),  ~,  ...
         S.Rm.(fn{i}){f,j},  ~  ] = glrt(A(:, idx.(lower(fn{i})){j}) ,C,K);
     
    end
end


end

function test_case
% DEV TEST
%
% ... and where is this from? 

tic
load /m_files/test_data/doa_on_range_cell.mat

CFG.dmth = {'mu'}; 
CFG.mus_param = [10 5 8]; % [20 10 3] or [40 20 2];
CFG.Nemit = n;

S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,CFG);

toc
keyboard

% test takes 0.010177 sec on optiplex 990

end




% DEPRECATED - DELETE LATER
% % this replaced by nan_pad (below) 
function err = store_err(idx,n)
% STORE DOA
% allows for possible return of 1 brg in 2 brg situation

% create storage
err = NaN(1,n);

% insert data
err(1:length(idx)) = idx; % sort(idx);


    
end

% this replaced by nan_pad (below) 
function doa_mus = nan_pad_doa(idx,th,n)
% STORE DOA - NaN pad DOA results
% allows for possible return of 1 brg in 2 brg situation
%
% 


% create storage
doa_mus = NaN(1,n);

% insert data ... used to be a sort here but I think the music solutions
% will get decoupled from the metrics if I do that
doa_mus(1:length(idx)) = th(idx); % sort(th(idx));


    
end
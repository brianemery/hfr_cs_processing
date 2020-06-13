function S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,n,dmth) %,snr)
% DOA ON RANGE CELL - run DOA methods on one rance cell of CS data
% S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,n,dmth)
%
% Typically called by doa_on_cs.m and experiment_music_vs_mle_roms.m - that
% is both by simulation code and by cs_processing code.
%
% *** run_param_test.m is disabled ***
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
% dmth    - (Optional) Cell struct of doa methods {'ml','mu','wf','sm','wm'}
%           defaults to {'ml','mu'}
%
% OUTPUTS
% Struct with radial sub-structs MU and ML based on MUSIC and MLE-AP, and more
% ... note that the bearing convention output is the same as that that is
% input (eg in APM) -- if it's cwN in, then it's cwN out. 
%
% SEE ALSO
% doa_on_cs.m, radial_from_cs.m, cs_processing.m, run_cs_processing*

% NOTE
% - check that error storage is correct? 

% TO DO
% I wonder if consolidating struct storage would be worth doing ...
% This function could be greatly simplified by the consistent use of one
% method to insert data into the DOA struct, eg using
% create_doa_column_index.m

% test?
if strcmp('--t',CS), test_case, return, end

% Enable ability to specify DOA method to use
if nargin < 7
    dmth = {'ml','mu'}; % default to just two
end


% get array matrix for music error (SeaSondes)
A = get_array_matrix(APM);

% get data storage structure
S = doa_struct(max([length(peakIdx) 1]),n,size(A,1),dmth);
S.Type = [APM.Type ' DOA Struct'];

mus_err = cell(n,1);



% loop over peak indicies
for f = 1:length(peakIdx) 
    
    % Get fbin for this iteration
    fbin = peakIdx(f);
    
    % build covariance matrix     
    C = make_cov(CS,fbin,rdx);
 
        
    % RUN DF METHODS
    % ... one est for each DOA solution ...
    [doa,ftime,idx,D,V] = run_doa_calculations(A,C,APM.BEAR,n,dmth); 
    

    % Get MUSIC Error Ests
    for i = 1:n      
        mus_err{i} = music_error(A,C,APM.BEAR,K,i,idx.mu{i});
        
        % pad with NaN if necessary
        mus_err{i} = store_err(mus_err{i},i);
        
    end
     
    
    % STORE OUTPUTS
    % 'ml' and 'mu' get converted to ML and MU here
    S = update_storage(S,doa,ftime,f,idx);
    
    
    % Save True Velocity
    S.RadVel(f,1:size(S.RadVel,2)) = CS.Vrad(fbin);
    
    % Save MUSIC things 
    S.Err(f,:) = [mus_err{:}];   
    S.eigValues(f,:) = D(:).'; 
    S.OtherMatrixVars.eigVectors{f,1} = V;
    
    % save this for investigations ...
    S.OtherMatrixVars.cov{f,1} = C;
    
    % seasonde music params
    if size(C,1) == 3
        [S.Params(f,:),S.Dual(f)] = run_param_test(D,V,A,idx.mu{2},[40 20 2]);
    end
    
    
    % compute likelihood ratios for all DOA methods
    S = calc_and_store_glrt(S,idx,A,C,K,n,f);
        
    % *** SIGNAL POWER HERE *** .
    % S = calc_and_store_power(S,idx,A,C,n,f);

    
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
S = calc_and_store_power(S,A,n);



end

function [doa,ftime,idx,D,V] = run_doa_calculations(A,C,th,n,fn)
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


if nargin < 5
    % fn = {'ml','mu','wf','sm','wm'};
    fn = {'ml','mu'}; 
end


% INIT OUTPUTS
for i = 1:numel(fn)
      idx(1).(fn{i}) = {};
    ftime(1).(fn{i}) = [];
    doa(1).(fn{i}) = {}; %NaN(1, (n*(n+1))/2 );
end


% COMPUTE DOAs, TRACK TIME

% MLE = CML = DML
if ismember('ml',fn)
    tic
    for i = 1:n
        idx(1).ml{i,1} = mle_ap(A,C,i);
    end
    ftime(1).ml = toc;
end


% MUSIC
if ismember('mu',fn)
    tic
    [~,idx.mu,D,V] = music(A,C,th,n);
    ftime.mu = toc;
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

for j = 1:numel(fn)
    for i = 1:n

        doa.(fn{j}){i,1} = nan_pad_doa(idx.(fn{j}){i},th,i);
        
    end
end

 
% reshape outputs
for i = 1:numel(fn)   
    doa.(fn{i}) =  [doa.(fn{i}){:}];    
end


end

function doa_mus = nan_pad_doa(idx,th,n)
% STORE DOA - NaN pad DOA results
% allows for possible return of 1 brg in 2 brg situation
%
% 


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

function S = calc_and_store_power(S,A,max_n) % idx,A,R,n,f)
% CALCULATE AND STORE POWER - deal with the DOA substructs
% S = calc_and_store_power(S,idx,A,C,n);
%
% this handles the substructs and getting the output in the right format
%
%
% S 	- Radial DOA structure
% idx   -
% A     - Antenna Array matrix (see array_matrix.m)
% R     - 
% max_n - Search for up to max_n emitters  
% f     -
%
% NOTE
% this function is designed to run on the whole struct, thus making it
% useful outside of this main function. Maybe this is a better way to do
% this kind of thing? rather than in the loop?

% get indexing to map the number of emitters to the column index
col = doa_column_index(max_n);

% get the doa method substruct field names
fn = fieldnames(S.Bear);

for n = 1:max_n                % number of emitters
    
    for r = 1:size(S.RadVel,1) % row index
        
        for i = 1:numel(fn)    % sub struct index
            
            % get the APM index of the DOA solutions
            ix = S.Idx.(fn{i}){r,n};
            
            if ~isempty(ix)
                
                % compute the power - this outputs an nxn matrix
                Pwr = signal_power( A(:,ix) , S.OtherMatrixVars.cov{r} );
                
                % get diagonal of the power matrix and put it in place
                %  Note that music might output less than n emiters so
                % we need to account for that possibility
                pwr = NaN(1,n);

                p = real(diag(Pwr));
                
                pwr(1:length(p)) = p;
                                
                S.Pwr.(fn{i})(r,col{n}) = pwr;
                
            end
            
        end
    end
end


end

function test_case
% DEV TEST
%
% ... and where is this from? 

tic
load /m_files/test_data/doa_on_range_cell.mat

S = doa_on_range_cell(CS,APM,peakIdx,rdx,K,n);

toc
keyboard

% test takes 0.010177 sec on optiplex 990

end


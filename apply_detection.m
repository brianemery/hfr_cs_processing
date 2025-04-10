function D = apply_detection(S,em,method,sfn) %,tf)
% APPLY DETECTION - apply detection result to extract DOA solutions
% S = apply_detection(S,em,method) %tf)
%
% Reform a DOA struct (containing Radial Data from a single time period)
% to keep only the DOA solutions that pass the given test (GLRT, SeaSonde,
% MUSIC-highest, etc.). 
%
% NOTE
% - this assumes that the RangeBearHead field is only one column containing
% the Range information. Bear is also a field and these will be merged
% later
% 
%
% 
% INPUT
% S  - DOA struct with LR field and subfields for DOA methods (eg. MU, ML)
% em - col vector giving the number of emitters in each row of the DOA
%      struct
% tf - boolean true if detection threshold was met (em might contain a
%      a value for the number of emitters even if threshold was not met)?
%
% method = char string for documenting the method used
%
% sfn  = option other non-standard field names to apply the single column
%        treatment to. 
%
% written for GLRT with the eye toward a general use case, eg any data with
% the same format (boolean) could be used to reshape the DOA data
%
% EXAMPLES
%     % compute emitters
%     [em,~] = compute_emitters_from_lr(LR,cut(i));
%     
%     % apply detection
%     D = apply_detection(S,em,'glrt'); 
%
%     % get SeaSonde detection
%      S = detection_codar_post_proc(S,APM);
% 
%     % run for seasonde
%     U = apply_detection(S,S.Dual+1,'seasonde'); 
%
% where cut appears to be 25 for the SS, and 300 for lera
% 
%
% Typically these run together:
% - get_radial_meta.m
% - apply_detection.m
% - doa_to_radial_struct.m or maybe detection_codar_post_proc.m
%
% SEE ALSO
% compute_emitters_from_lr.m, detection_codar_post_proc.m, run_param_test.m
%
% get_radial_meta.m ususally runs before this, see rng_processing_bs.m's
% subfunction cs_processing, ... then it's followed by
% doa_to_radial_struct.m ... jeesh

% TO DO
% ... need to compare these and make sure they get the same thing ...
% eg: D = apply_detection(S,S.Dual+1);
%  and apply_test_result.m which works for the M=3 case (e.g. SeaSonde)
% 
% detection_codar_post_proc.m ...
%
% NEED TO PUT THE LR, or other detection info in OtherMetadata struct

% TO DO
% want to be able to vary the LR threshold ... 
%
% Need to be able to have within the struct different detection results, eg
% this applies only to MU or ML, etc ...

% check if we're running the test
if strcmp(S,'--t'), test_case, return, end


% add emitter field to the input struct to propagate it through
S.Nem = em;

% get output struct to modify, also eject if no emitters 
D = S; if isempty(em), return, end

% create a cell array that maps the number of emitters to the columns of
% the DOA matricies
cols = doa_column_index(5);




% Create a boolean array to use to extract the data we want to keep based
% on the determined number of emitters.

% allocate the boolean
tf = false(size(S.RadVel));

% populate the ones to keep with 'true'. This uses the number of emitters
% in 'em' to get the column index. Handle the odd case where we have NaN
% emitters (or zero). This passes them through as all false, which gets
% handled appropriately by the logical indexing below.
for i = 1:size(tf,1)
   
    if ~isnan(em(i)) & em(i)~=0
        tf(i,cols{em(i)}) = true;
    end
    
end



% Apply the boolean to struct field elements 

% apply to top level fields 
fn = {'RadVel','Err','Apprch','PkPwr','PkWdth','eigValues'}; %,'Params'};
fn = fn(ismember(fn,fieldnames(S)));

% note that the application of tf here creats a column vector in D, unless
% the input is a vector!
for i = 1:numel(fn)
    D.(fn{i}) = S.(fn{i})(tf);
end

% struct fields that contain DOA substructs
fn = {'Bear','RmsBrgErr','BrgDiff','RomsBrg','Pwr','Pwr2','Pwr3'}; 
fn = fn(ismember(fn,fieldnames(S)));

for i = 1:numel(fn)
    D.(fn{i}) = apply_to_substruct(S.(fn{i}),tf);
end



% These get the single column treatment
%
% This arbitrarily has taken the SNR of element 3 - for the ULA and RA8
% there might be a better way ... 
fn = {'SNR','Dual','RngIdx','Nem','RangeBearHead'}; 
if nargin > 3, fn = [fn sfn]; end 
fn = fn(ismember(fn,fieldnames(S))); 

D = single_col_treatment(S,D,fn,tf);

if isfield(S,'OtherMatrixVars')
    fn = {'eigVectors','cov'};
    fn = fn(ismember(fn,fieldnames(S)));
    
    D.OtherMatrixVars = single_col_treatment(S.OtherMatrixVars,D.OtherMatrixVars,fn,tf);
end

% Special treatment for Params
% Take each column, repmat it into a 3 column matrix, apply tf, then
% reconstruct the original matrix format
if ~isempty(S.Params)
    dm1 = repmat(S.Params(:,1),1,size(tf,2));
    dm2 = repmat(S.Params(:,2),1,size(tf,2));
    dm3 = repmat(S.Params(:,3),1,size(tf,2));
    
    D.Params = [dm1(tf) dm2(tf) dm3(tf)];
end


% SPECIAL CASE HANDLING
% Deal with the logical indexing 'bug', which is only encountered if S
% contains results from a single Doppler bin
%
% ... maybe just define a col function for everything ...
if size(S.RadVel,1) == 1
    
    
    % for this, define a colon function to force column vectors
    % see note about logical indexing bug next comment block
    col = @(x) x(:);
    D.Params = [col(dm1(tf)) col(dm2(tf)) col(dm3(tf))];
    
    
    % all the fields need to be turned into column vectors in this case
    fn = {'RadVel','Err','SNR','Dual','RngIdx','RangeBearHead','Nem'};
    
    for i = 1:numel(fn)
        D.(fn{i}) =  D.(fn{i})(:);
    end
    
    % struct fields that contain DOA substructs
    fn = {'Bear','RmsBrgErr','BrgDiff','RomsBrg'};
    fn = fn(ismember(fn,fieldnames(S)));
    
    for i = 1:numel(fn)
        D.(fn{i}) = apply_to_substruct(D.(fn{i}),':');
    end
        
end

% NOTE ON BELOW
% I think the solution is the expand the LR fields to the columns in col
% (same as Bear arrays), and then treat same as those ones


% -- STATUS --
% REGARDING ATTEMPT TO KEEP LR VALUES ...
% I'm going to consider this for now as a deep dive into the weeds, that is
% tangent to the main focus of the paper. Tangent, but important, but let's
% make sure the main focus gets taken care of first BE 4/4/2019
% 
% 
% % KEEP LR VALUES
% % might need to check this more carefully since this is being done at
% % about 50% mental capacity ...
% %
% % ... I think what I am doing is getting the vector eg from the LR matrix
% % that only contains the LR values that were used ... then expanding this
% % into a matrix that can be used with tf (in the same way the single col
% % data is above), then applying tf ...
% %
% % these are #rows by cols= emitter number ... seems to work for cell also
% fn = {'LR','GM','LR2','Rm'};
% 
% % REDUCE 
% % This makes the single column equivalent
% 
% % convert the row and column index implied by em into an absolute index
% ix = sub2ind(size(S.LR.ML),(1:length(em))',em);
% 
% for i = 1:numel(fn)
%     D.(fn{i}) = apply_to_substruct(S.(fn{i}),ix);
% end
% 
% % Save number of emitters ... use single column treatment method, sort of
% Em = repmat(em,1,size(S.LR.ML,2));
% 
% D.Em = Em(ix);
% 
% % Now expand this to column size that the tf can be applied to ...
% 
% % dummy variable, created from the D struct data from above
% dm = repmat(D.Em,1,size(tf,2));
% 
% D.Em = dm(tf);
% 
% % Do this for the substructs ...
% keyboard
% fn = {'LR','GM','LR2','Rm'};
% 
% % for i = 1:numel(fn)
% %     
% %     % dummy variable, created from the D struct data from above
% %     dm = repmat(D.(fn{i})(:,1),1,size(tf,2));
% %     
% %     D.(fn{i}) = dm(tf);
% %     
% % end
% 


% will have to sort these out later ...
% fn = {'eigValues'} 


% document
if nargin > 2, D.ProcessingSteps{end+1} = method; end

% include this mfile
D.ProcessingSteps{end+1} = 'apply_detection';

% idempotent cleanup
D.LonLatUnits = {'Decimal Degrees','Decimal Degrees'};
D.RangeBearHeadUnits = {'km','Degrees_ccw_from_east','Degrees_ccw_from_east'};




end

function D = apply_to_substruct(S,tf)
% APPLY_TO_SUBSTRUCT
%
% apply detection to DOA method substructs

D = S;

if isstruct(S)
        
    fn = fieldnames(S);
    
    for i = 1:numel(fn)
        if ~isempty(S.(fn{i}))
            D.(fn{i}) = S.(fn{i})(tf);
        end
    end

end

end

function D = single_col_treatment(S,D,fn,tf)
% fn = {'SNR','Dual','RngIdx','Nem','RangeBearHead'}; 

for i = 1:numel(fn)
    
    dm = repmat(S.(fn{i})(:,1),1,size(tf,2));
    
    D.(fn{i}) = dm(tf);
    
end

end

function test_case
%
% TEST CASE
%
% discovered a 'bug': application of logical indexing will create a column
% vector when it's a matrix and the variable is a matrix. However, it will
% create a row vector when the variable is a row vector. In general,
% logical indexing *ALWAYS* turns matricies into column vectors, and
% produces a vector the same orientation as it starts with. 

% DATA FROM running run_cs_post_processing.m on
% some of Anthony's LERA data, where the logical indexing
% bug was encountered. 
% 
% normal data set to test general functionality
load /m_files/test_data/apply_detection_1.mat

D = apply_detection(S,em,method);

keyboard
 
% this one has a single range cell, single point with a two emitter
% solution
load /m_files/test_data/apply_detection_2.mat

D = apply_detection(S,em,method);

keyboard

% % logical indexing test code
% m = reshape(1:9,3,3);
% tf = false(size(m));
% tf([1 5 7]) = true;
% v = 1:9
% 
% v(tf)
% 
% v = v(:)
% 
% v(tf)
% 
% tf = tf(:);
% m(tf)
% 
% tf = true(size(m))
% 
% m(tf)


end
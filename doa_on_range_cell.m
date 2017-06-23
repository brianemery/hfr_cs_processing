function [MU,ML] = doa_on_range_cell(CS,APM,peakIdx,rdx,K) %,snr)
% DOA ON RANGE CELL - quick and dirty comparison
% [MU,ML] = doa_on_range_cell(CS,APM,peakIdx,rdx)
%
% Typically called by doa_on_cs. for example ...
% Make sure CS have not been converted to dBm
%
% INPUTS
% ... K snapshots for music error
%
% OUTPUTS
% radial structs MU and ML based on MUSIC and MLE-AP
%
% SEE ALSO
% doa_on_cs.m, radial_from_cs.m, cs_processing.m
%
% NOTE
% Gets used by some of the radar simulation code so it's still around. For
% the processing of real Cross Spectra, see radial_from_cs.m or similar


% INIT OUTPUTS
% matrix of single and dual bearing solutions for both methods
% index refers to the index of the APM
ML = doa_struct(length(peakIdx)); 
MU = ML;

MU.Type = 'MUSIC';
ML.Type = 'MLE-AP';

% get array matrix for music error (SeaSondes)
A = get_array_matrix(APM);


% loop over peak indicies
for f = 1:length(peakIdx);
    
    fbin = peakIdx(f);
    
    % build covariance matrix     
    C = make_cov(CS,fbin,rdx);
    
 
    % RUN DF METHODS
    
    % Max liklihood - alt projection NEED TO CHECK CODE - rm ends?
    [ix,~] = mle_ap(APM,C);

    % Max liklihood
    % [ix,~] = mle(APM,C);
    
    % Music
    % keyboard
    % [~,sdx,ddx,tr,P] = music(APM,C);
    [~,mx,D,V] = music(A,C,APM.BEAR);
    
    % .. could use a clean up
    [P,tr] = run_param_test(D,V,A,mx{2});
    
    % stack indecies in a row vector
    % account for possibility of length(dualIdx) = 1
    mx = [mx{1} mx{2}']; %[sdx ddx(:)'];
    
    % music
    MU.Bear(f,1:length(mx))   = APM.BEAR(mx);
    MU.RadVel(f,1:length(mx)) = CS.Vrad(fbin);
 
    % mle
    ML.Bear(f,:)   = APM.BEAR(ix);
    ML.RadVel(f,1:length(ix)) = CS.Vrad(fbin);
    
    
    % Just use MUSIC Parameter test for now ...
    MU.Dual(f) = tr;
    MU.Params(f,:) = P;
    ML.Dual(f) = tr;
        
    % Guess these should take V,D as inputs ...
    MU.mus_err(f,1) = music_error(A,C,APM.BEAR,K,1,mx(1));
    MU.mus_err(f,2:length(mx)) = music_error(A,C,APM.BEAR,K,2,mx(2:end));
    
%     MU.mus_err2(f,1) = music_error(A,C,APM.BEAR,K,1,sdx,snr(f,1));
%     MU.mus_err2(f,2:length(ddx)+1) = music_error(A,C,APM.BEAR,K,2,ddx,snr(f,1));
%     
    %     % PLOTTING
    %
    %     % MUSIC dual and single
    %     hu    = plot(MU.Bear(f,2:3),MU.RadVel(f,2:3),'go');
    %     hu(2) = plot(MU.Bear(f,1),MU.RadVel(f,1),'g*');
    %
    %     % MLE dual and single
    %     hx    = plot(ML.Bear(f,2:3),ML.RadVel(f,2:3),'ro');
    %     hx(2) = plot(ML.Bear(f,1),ML.RadVel(f,1),'r*');
    %
    %
    %     if f == 1,
    %         legend([hu(:)' hx(:)'],'MUSIC dual','MUSIC sngl','MLE dual','MLE sngl')
    %
    %     end
    %
    %     pause
end







end





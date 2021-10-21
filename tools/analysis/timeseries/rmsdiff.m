function [rmsd,bias,N,mae]=rmsdiff(a,b)
% RMSDIFF.M - RMS differences 
% [rmsd, bias, N, mae]=rmsdiff(x,y);
% Computes the unbiased root-mean-square difference
% between 2 timeseries. This function assumes that 
% input data time increments columnwise, allowing it
% to accept matricies. Bias's, the difference of the
% means, is output.
%
% Also outputs the old rms diff for comparison, and now
% and the mean absolute error (MAE)
%
% see also rsquared.m, mean_noNaN.m
%
% CONFIDENCE INTERVALS
% 
%     % rmsd ci? ... using matlab's bootstrap 
%     [ci,bootstat]  = bootci(1000,@rmsdiff,x,y);
%     [ci,bootstat]  = bootci(1000,{@rmsdiff,x',y'},'Type','norm');
%     
%     rmsd_ci(1,i) = ci(2);
%     rmsd_ci(3,i) = ci(1);
%     rmsd_ci(2,i) = rmsdiff(x,y);


% 20mar99 Brian Emery
% 21Feb01 Updated: Remove means, report bias's,
% Matrixified the code.

% get rid of NaN's. ad and bd have time increment row-wise
ad=a'; bd=b';
c=ad+bd; 
i=isnan(c); 
ad(i)=0; bd(i)=0;
abar=(sum(ad,1)./sum(~i,1))'; bbar=(sum(bd,1)./sum(~i,1))';

% compute N
N=sum(~i,1);

% compute bias
bias=abar-bbar;

% Re-insert NaN's, remove means and compute the difference of the 2 timeseries's
ad(i)=NaN; bd(i)=NaN;
del=(ad-(abar*ones(1,size(ad,1)))')-(bd-(bbar*ones(1,size(bd,1)))');  

% compute rms diff, again using careful treatment of NaN's
del(i)=0;
rmsd=sqrt(sum(del.^2,1)./N)'; % previously ./(sum(~i,1))

% compute mae https://en.wikipedia.org/wiki/Mean_absolute_error
mae = (sum(abs(del),1)./N)';



% % output rms from old method for comparison
% rmsold=oldmethod(a,b);

return


%%%%
% alternative method of computing rmsDiff:
rmsd2=[]; rmsdiff=[];

for i=1:row
	% get rid of NaN's, j is the column index
	j=find(~isnan(a(i,:)+b(i,:)));

	del=a(i,j)-b(i,j);

	[xbar,stdev,hi,lo]=stats_noNaN(del);
   rmsdiff=sqrt((stdev.^2)+ (xbar.^2));
   rmsd2=[rmsd2 rmsdiff];
end
end


function rmsold=oldmethod(a,b)
% OLD METHOD: 
% run a loop one row at a time
[row,col]=size(a);
rmsd=[];

for i=1:row
	% get rid of NaN's, j is the column index
	j=find(~isnan(a(i,:)+b(i,:)));

	del=a(i,j)-b(i,j);

	rmsdiff=sqrt( (1/(length(del)-1))*sum(del.^2)  );

	rmsd=[rmsd rmsdiff];
end

rmsold=rmsd;
end

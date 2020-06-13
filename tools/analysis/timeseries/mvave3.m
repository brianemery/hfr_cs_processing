 function flt = mvave3(dat,wn)
% MVAVE3.M
%  [flt]=mvave3(data,i);
%  Filter data using a moving average. Where:         
%        -'data' is the data to be filtered (vector or matrix)
%        -'i' is a number which gives the size of 
%         the filter (7 for tow 6, and tow 2 salinities)
%        -flt is the filtered data 
%
% This function computes a moving average for each row if 'data'
% is a matrix

% Copyright 1995-2010 Brian Emery
%
% Version 1 Written by Brian Emery, sometime in 1995, using Libe 
% Washburn's idea.


flt = dat; % this allows you to keep the ends

% half window length
wi = round(wn/2);

% loop over rows
for k = 1:size(dat,1)
    
    flt(k,wi:end-wi) = row_mv_ave(dat(k,:),wn);
    
%     % name the current row 'data'
%     data = data2(k,:);
%     
%     % initthe matrix to be averaged
%     D = zeros(i,c+(i-1));
%     
%     % populate the matrix to be averaged
%     for m = 1:i
%         D(m,m:c+(m-1))=data;
%     end
%     
%     % filter the data using mean to take the average
%     flt = mean_noNaN(D(:,i:c).').';
%     
%     % build the matrix back up
%     filter2(k,:) = data;
%     filter2(k,  ) = flt; keyboard
      
end

 end

 function flt = row_mv_ave(dat,wn)
 % Compute Moving Average of just one row vector
 
 c = size(dat,2);
 
 % initthe matrix to be averaged
 D = zeros(wn,c+(wn-1));
    
 % populate the matrix to be averaged
 for m = 1:wn
     D(m,m:c+(m-1)) = dat;
 end
 
 flt = mean_noNaN(D.').';
 
 % only keep meaningful stuff
 flt = flt( (wn-1):(end-wn-1));
 
 
 end
 
 
 
 
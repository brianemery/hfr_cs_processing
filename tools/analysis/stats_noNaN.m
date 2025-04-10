function [xbar,stdev,hi,lo,med,n]=stats_noNaN(x)
% STATS_NONAN
% [xbar,stdev,hi,lo,median,n] = stats_noNaN(x)
% Computes the mean, std, and ranges of each row of a
% matrix after removing NaN's. Also returns number of
% non-NaN values used in the calc.
%
% SEE ALSO
% mean_noNaN.m
% make_stats_table.m

% Copyright (C) 1999-2010 Brian M. Emery
% Brian Emery 13Sept99
% try-catch added 3apr09

% check for test case
if strcmp('--t',x), test_case, return, end


[r,c]=size(x);

% if x is row vector, orient it so that it computes the
% mean of the vector
if r==1 || c==1
   x=x(:)';
   r=size(x,1);
end

% create empty outputs
[xbar,stdev,hi,lo,med,n]=deal(NaN(r,1));

% If there is a way to do this without a loop, I'd like to know what it is.
for i=1:r
        row = x(i,~isnan(x(i,:)));
        
       n(i) = length(row);
    xbar(i) = mean(row);
   stdev(i) = std(row);
     med(i) = median(row);

  % when NaN's are put in, these return empty matricies and an error
  try hi(i) = max(row); 
  catch
      hi(i)=NaN; 
  end
  
  try lo(i) = min(row); 
  catch
      lo(i)=NaN; 
  end

end



end

function test_case
% test cases:
%       Generate normal values with mean 1 and standard deviation 2.
%
n = 10000; 
x = 1 + 2.*randn(n,1);

plot(x,'-b.')
hold on


% Generate 10% NaNs, use random integers uniform on the set 1:n 
n = ceil(100.*rand(n/10,1));
 

% insert some NaN's
x(n)=NaN;

plot(x,'r.')

tic
[xbar,stdev,hi,lo,med]=stats_noNaN(x);
toc

keyboard


return


clear
% make x a matrix to test that too
x(1,:)= 1 + 2.*randn(100,1);
x(2,:)= 2 + 4.*randn(100,1);
x(3,:)= 3 + 6.*randn(100,1);
x(4,:)= 4 + 8.*randn(100,1);
x(5,:)= 5 + 10.*randn(100,1);

n = ceil(500.*rand(50,1));
x(n)=NaN;

end	

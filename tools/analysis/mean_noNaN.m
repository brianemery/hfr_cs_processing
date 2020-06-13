function abar=mean_noNaN(a)

% MEAN_NONAN
% abar=mean_noNaN(a)
% Computes the mean of each row of a
% matrix after removing NaN's. 

% see also mean_noNaN_old.m

% Copyright (C) 1999-2010 Brian Emery
% Updated 5Oct1999  

ad=a';
i=isnan(ad);
ad(i)=0;
abar=(sum(ad)./sum(~i))';
end



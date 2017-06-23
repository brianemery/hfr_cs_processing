function T = hourly_merge(R)
% HOURLY MERGE - paste together other mfiles to merge 10 min data to hrly
% Create info needed for binning
%
% ... need to look at this a bit more carefully ...


% bin width (min convert to days)
bw = 72/1440;

% bin centers
bc = min( round(R.TimeStamp*24)/24 ):(1/24):max( round(R.TimeStamp*24)/24 );


T = bin_data_struct(R,'TimeStamp',bw,bc);



end

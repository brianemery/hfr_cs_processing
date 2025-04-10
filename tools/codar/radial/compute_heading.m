function R = compute_heading(R)
% COMPUTE HEADING - from site origin and data lon lat using dist.m
%
% The dist.m output is degrees cwN, which this converts to ccwE for use
% with HFR Progs
% 
% EXAMPLE
% R = compute_heading(R)

% Should combine this with rangeBear2LonLat.m

% Brian Emery, 5 Aug 2019

% check for test 
if strcmp('--t',R), test_case, return, end


% struct recursion if necessary
if numel(R) > 1
    for i = 1:numel(R)
        R(i) = compute_heading(R(i));
    end
    return
end

% Run the code 
[R.RangeBearHead(:,3),af]= calc_head(R.SiteOrigin,R.LonLat);


% convert to ccwE
R.RangeBearHead(:,3) = cwN2ccwE(R.RangeBearHead(:,3));

% % set units %THIS DOESNT WORK WITH the recursion code
% if ~isfield(R,'RangeBearHeadUnits')
%     R.RangeBearHeadUnits ={'km','?','deg ccwE'};
% else
%     R.RangeBearHeadUnits{3} = 'deg ccwE';
% end


end

function [Head,af] = calc_head(SiteOrigin,LonLat)
% [Head,ar] = calc_head(SiteOrigin,LonLat)
%
%       site to grid direction = bearing_f = Bear
%       grid to site direction = bearing_r = Head



[Head,af] = deal(NaN(size(LonLat,1),1));

SiteOrigin = SiteOrigin(1,:);


% just loop I guess
for i = 1:size(LonLat,1)
    
    % [RANGE,AF,AR]=dist(LAT,LONG)
    % Site origin is Lon Lat
    [~,af(i),Head(i)] = dist([SiteOrigin(2) LonLat(i,2)],[SiteOrigin(1) LonLat(i,1)]);
    
    
end
    
    
end

function test_case
% TEST CASE
% 
% ... tested during the development of submesoscale_simulation




end

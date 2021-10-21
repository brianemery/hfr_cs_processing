function cwN=ccwE2cwN(ccwE)

% CCWE2CWN.M
% cwN=ccwE2cwN(ccwE)
% Converts Degress Counter-clockwise from east to 
% clockwise from north (ie, from cartesian to true bearing)

% 15mar99 Brian Emery

cwN=450-ccwE;
i=find(cwN>=360);
cwN(i)=cwN(i)-360;


% NOTE from wavemeas.pdf
% in radians:
% cwn = 3pi/2 - ccwE


end


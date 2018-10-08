function ccwE=cwN2ccwE(cwN)

% CWN2CCWE.M
% ccwE=cwN2ccwE(cwN)
% Converts Degrees clockwise from North to Degrees
% counter-clockwise from East (CCW E). 

% Brian Emery 1Feb98 

% This is the ugly, brute force method.
ccwE=(abs(cwN-360))+90;
i=find(ccwE>360);
ccwE(i)=ccwE(i)-360;
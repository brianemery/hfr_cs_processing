function deg = rad2_to_deg(rads2)
% RAD^2 TO DEG - radians squared to degrees
% ... converts the variance to a standard deviation
%
% c.f. music_error.m, music_error2.m

deg = sqrt( rads2.* ((180/pi)^2) );


end
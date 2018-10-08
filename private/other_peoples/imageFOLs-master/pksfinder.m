function [pks,dzdt] = pksfinder(db,db_level);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [pks,dzdt] = pksfinder(db,db_level);
%
%  same as updown.m, renamed for HFR_DP
%
% Creates an indexing array marking when a moving profiler has changed
% directions.  Also creates a flag array, equal in length to the
% original data, indicating the direction of travel.  Note, the direction
% is relative to the pressure record's syntax.
%
% Inputs:
%	db		--	depth/pressure record
%	db_level	--	threshold indicating a "true"
%				change in direction.
%
% Outputs:
%	pks		--	indexing array marking when the direction of
%				travel (up or down) changes.
%	dzdt		--	flag array indicating the direction the
%				profiler is moving. equal in length to db.
%				note: dzdt is relative to pressure record
%				syntax.
%
% Example:
%	[pks,dzdt] = updown(ctd(:,2),10);
%
% Originally written by:
%	Christopher Wingard
%	Sr. Faculty Research Assistant
%	College of Oceanic and Atmospheric Sciences
%	Oregon State University
%	Corvallis, Oregon
%	March 2003
%
% Editted by (slower, but more accurate):
%	Anthony Kirincich
%	Ph.D. Candidate
%	College of Oceanic and Atmospheric Sciences
%	Oregon State University
%	Corvallis, Oregon
%	January 2005
%
% Re-Editted by Chris Wingard (faster and more accurate) March 2005
%
% Special note:
%	Calls the MM6 Toolbox function mmpeaks.
%	For more info see http://www.eece.maine.edu/mm/.
%	-- I coppied the relevant code into here (small violation I know, don't
%	-- shoot me!) so others could use this without having to get the mm6 
%	-- toolbox. 

% determine the up/down cycle of the profiler and mark the local min/max extrema
% via an indexing array (code from mmpeaks.m).
% pks = mmpeaks(db,'all');	% note pks equals an indexing array
y = db; ys = size(y); x = reshape(y,1,prod(ys));
j = sign(diff([-inf x -inf]));
j = find(diff(j+(j==0)) == -2);
k = sign(diff([inf x inf]));
k = find(diff(k+(k==0)) == 2);
pks = sort([k j])'; clear y ys x j k

%%%%%%%%% new code as of March 2005 from cwingard %%%%%%%%%

% mark and delete index locations where there isn't a "true" direction change
cnt = 1; step = 1; m = [];
while cnt < length(pks)
	if cnt + step == length(pks)
		break
	end %if
	% compare consequtive extrema
	dz1 = db(pks(cnt)) - db(pks(cnt+step)); % magnitude of difference
	sg1 = sign(dz1);
	if abs(dz1) < db_level % change is not large enough
		m = [m; cnt + step]; % mark extrema for deletion
		step = step + 1; % bump the stepper
	else % change is large enough, but does it qualify as a direction change?
		if cnt + step == length(pks)
			break
		end %if
		% compare the next set of extrema to see if this is a real change
		dz2 = db(pks(cnt+step)) - db(pks(cnt+step+1)); % magnitude
		sg2 = sign(dz2); j = cnt + step;
		if abs(dz2) < db_level % change is not large enough -- or is it ???
			if sg1 == sg2 % nahh, we're still going the same way
				m = [m; cnt + step]; % mark extrema for deletion
				step = step + 1; % bump the stepper
			else % maybe -- this slows things down
				% find extrema > db_level to the left of the current point
				lft = db(pks(1:j-1)) - db(pks(j));
				i = find(abs(lft) > db_level);
				if isempty(i) == 1
					a = 1;
				else
					a = i(end);
				end
				% find extrema > db_level to the right of the current point
				rght = db(pks(j+1:end)) - db(pks(j));
				i = find(abs(rght) > db_level);
				if isempty(i) == 1
					b = length(pks);
				else
					b = j + i(1);
				end
				% set the min and max of the window
				mn = min(db(pks(a:b)));
				mx = max(db(pks(a:b)));
				% see if our point is a min or max of this window
				if mx == db(pks(j)) | mn == db(pks(j))
					% we have a valid turn -- reset the counters
					cnt = j; step = 1;
				else % it's not a valid peak
					m = [m; j]; % mark extrema for deletion
					step = step + 1; % bump the stepper
				end %if
			end %if
		else % we have a valid turn -- reset the counters
			cnt = j; step = 1;
		end %if
	end %if
end %while
pks(m) = []; clear cnt step j lft rght i a b mn mx m dz

% mark the tail
if pks(end) ~= length(db)
	pks = [pks; length(db)];
end %if

flg = ones(length(pks),1);
for i = 2:length(pks) - 1
	flg(i) = (db(pks(i)) > db(pks(i-1)) & db(pks(i)) > db(pks(i+1))) | ...
		(db(pks(i)) < db(pks(i-1)) & db(pks(i)) < db(pks(i+1)));
end
pks(flg==0) = []; clear flg i

% create the direction array
dz = sign(diff(db(pks))); dzdt = zeros(length(db),1);
for i = 1:length(pks)-1
	dzdt(pks(i):pks(i+1)) = dz(i);
end
clear dz i

return
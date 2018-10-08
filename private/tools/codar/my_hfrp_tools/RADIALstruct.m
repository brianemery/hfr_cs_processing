function R = RADIALstruct( N )
% RADIALSTRUCT  Creates an empty default RADIAL structure 
% that can be filled with radial current measurements data
%
% Usage: RADIAL = RADIALstruct( N )
%
% This function should be used *every* time one wants to load in any type of
% radial data so that all such structures have the same fields.
%
% Inputs
% ------
%
% N: size of vector of radial structures to return.  Defaults to 1.
%
% Outputs
% -------
%
% RADIAL: is the empty structure to use for recording radial data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: RADIALstruct.m 396 2007-04-02 16:56:29Z mcook $	
%
% Copyright (C) 2007 David M. Kaplan
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BE Edit:
% added specific fields to the 'OtherMatrixVars' substruct

R.Type = 'Ideal';

R.SiteName = '';
R.SiteCode = 1;
R.SiteOrigin = [NaN,NaN];

R.FileName = cell(1,0); % This is the only way to have a truly empty RADIAL
                        % structure
R.TimeStamp = zeros(1,0);
R.TimeZone = 'GMT';

R.LonLat = zeros(0,2);
R.RangeBearHead = zeros(0,3);

R.RadComp = [];
R.Error = [];
R.Flag = [];
R.U = []; % These are purely for convenience, not calculations
R.V = []; % These are purely for convenience, not calculations

R.LonLatUnits = {'Decimal Degrees','Decimal Degrees'};
R.RangeBearHeadUnits = {'km','Degrees_ccw_from_east','Degrees_ccw_from_east'};
[R.RadCompUnits,R.ErrorUnits,R.UUnits,R.VUnits] = deal( 'cm/s' );

R.CreateTimeStamp = datestr(now);
R.CreateTimeZone = 'GMT';


% Note that all variables in this structure should have the
% same number of dimensions as RadComp (i.e. [SpatialPoints x TimeSteps] ).
[R.OtherMatrixVars.ERSC , ...
 R.OtherMatrixVars.ERTC , ...
 R.OtherMatrixVars.ESPC , ...
 R.OtherMatrixVars.MAXV , ...
 R.OtherMatrixVars.MINV , ...
 R.OtherMatrixVars.SPRC , ...
 R.OtherMatrixVars.VFLG ]= deal([]); 


R.OtherMetadata = struct([]); 

R.ProcessingSteps = {};

R.RADIAL_struct_version = which(mfilename); %'SVN $Rev: 397 $ $Date: 2007-04-02 $';

if exist('N','var')
  R = repmat( R, [N,1] );
  
  % BE note: if you have N sites, this applies a site code 2^N for each
  % site. Ive been useing this function such that N may equal the number of
  % radial files I'm loading.
  % n = num2cell( 2.^(0:N-1) );
  % [R.SiteCode] = deal(n{:});
end

  
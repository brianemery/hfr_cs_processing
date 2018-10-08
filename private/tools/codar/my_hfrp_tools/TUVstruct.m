function T = TUVstruct( DIM, varargin )
% TUVSTRUCT  Creates an empty default TUV structure that can be
% filled with any type of total currents data
%
% Usage: TUV = TUVstruct( DIM, N )
%
% This function should be used *every* time one wants to generate any
% type of total currents (normal, OMA, etc.) as this will guarantee that
% they have the same structure.
%
% Inputs
% ------
%
% DIM: Size of U, V, etc. matrices to initialize.  A two element vector
% with [ NGridPts, NTimeStamps ].  Defaults to [ 0 0 ].
%
% N: The number of TUVerror structures to create in ErrorEstimates.  See
% TUVerrorstruct for details.
%
% Outputs
% -------
%
% TUV: is the empty structure to use for recording totals data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: TUVstruct.m 396 2007-04-02 16:56:29Z mcook $	
%
% Copyright (C) 2007 David M. Kaplan
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist( 'DIM' , 'var' )
  DIM = [0 0];
end

% Basics
T.Type = 'TUV';
T.DomainName = 'UCSB-SLO';

% BE Additions
T.SiteNames = {};
T.SiteOrigins = [];

T.CreationInfo = creation_info;  % BE Mod

% Time
T.TimeStamp = NaN(1,DIM(2));
T.TimeZone = 'GMT';

T.CreateTimeStamp = now;
T.CreateTimeZone = 'GMT';

% Space
T.LonLat = NaN(1,DIM(2));
T.Depth = NaN(DIM(1),1 );

% UV, linear combo of site codes, and number of radials (BE added)
[T.U,T.V,T.SiteUsed,T.NumRads] = deal( NaN(DIM) );

% boolean nSites x NTimeStamps? Not used yet. Need to see how HFRP was
% planning to do this
T.SiteCodes = [];


% Some units
T.LonLatUnits = { 'Decimal Degrees', 'Decimal Degrees' };
[T.UUnits,T.VUnits] = deal( 'cm/s' );
T.DepthUnits = 'm';

% Errors
T.ErrorEstimates = TUVerrorstruct( DIM, varargin{:} );


% Other
T.OtherMatrixVars = [];   % Should eventually be a structure.
T.OtherSpatialVars = [];  % Should eventually be a structure.
T.OtherTemporalVars = []; % Should eventually be a structure.
T.OtherMetadata = [];     % Should eventually be a structure.

T.ProcessingSteps = {};
T.OtherMetadata.TUV_struct_version = '/m_files/tools/codar/my_hfrp_tools/TUVstruct.m'; % BE Mod

end

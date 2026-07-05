%======================================================================
%                    M Y C R U I S E . M 
%                    doc: Wed Jun  4 12:22:17 2008
%                    dlm: Wed Jun  4 12:42:33 2008
%                    (c) 2008 M. Visbeck & A.M. Thurnherr
%                    uE-Info: 27 58 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

function [] = cruise(directory, dods)
% function [] = cruise(directory, dods)
%
% Create Cruise netCDF file from cast files.
%
% Input :  directory : Either a directory '/home/someone/lookhere/'
%                      Or the 'f' structure from LADCP where
%                      'f.res' will be used.
%
%          dods      : Structure with DODS info ** OPTIONAL **
%                      dods.cruise - dods address for cruise file
%                      dods.d1     - start of address for cast
%                      dods.d2     - end of address for cast
%
% Output : cruise.nc : Cruise netCDF file (filename will vary)
%          htmltable is called to create tables.

% Modifications by ANT:
%	Jun  4, 2008: - replace all spaces by _ in cruise name

if nargin < 1, help cruise, return, end
if nargin < 2, dods = ' '; end

% Interpolation Grids for Cruise netCDF File
GRID.z       = [10:10:6000]';      % Depth (m)
GRID.z_sadcp = [0:10:800]';        % Depth (m)
GRID.zbot    = [0:10:300]';        % Height Above Bottom (m)
GRID.tim     = linspace(0,1,100)'; % Normalized Time

% Attributes will be grouped together by first four letters:
AttVar = {'GEN_'; 'BAR_'; 'INPU'; 'LADC'};

% Main Varables, all others will be grouped by variable length as 'EXTRA'
MainVariables = {'z';'u';'v';'uerr';'lat';'lon';'jul';'maxz'; ...
                 'cast'; 'cast_id'; 'station'};

disp('CRUISE')
disp(' Surveying Files')
warning off

if isstruct(directory)
   loc = directory.res;
   ii = max(find(loc == '/'));
   loc = loc(1:ii);
else
   loc = directory;
end
d = dir([loc '/*.nc']);

count = 0;
bad = [];
for i = 1:length(d)
   f = netcdf([loc '/' d(i).name]);
   % Is this an LADCP File?
   if isempty(f.GEN_Profile_start_decimal_day(:))
      bad = [bad i];
   else
     count = count + 1;
   end
   close(f)
end
d(bad) = [];

%%%
%disp('TESTING - only using first 4 files')
%pause(2)
%d(5:end) = [];
%count = min(count, 4);
%%%

if count
   disp([' Found ' int2str(count) ' Files to Process']);
%   for i = 1:length(d)
%      disp(['  ' d(i).name]);
%   end
else
   error([' Nothing To Do!';])
end
ncquiet
f = netcdf([loc '/' d(1).name]);
name = regexprep(f.cruise_id(:),' ','_');
if isempty(name)
   name = 'unknown';
end
close(f);
disp([' Cruise ID: ' name ]);

% A : Attributes
% D : Data
% F : Dimensions and lists of Variables
% G : Dimension Lengths 
% H : Variable Type (char, float, etc.)

disp([' Surveying Files For All Possible Dimensions and Variables']);
[F,G,H] = deal([]);
for i = 1:length(d)
   d(i).nc    = d(i).name;
%   disp(['  ' d(i).nc ]);
   d(i).name  = d(i).name(1:end-3); % .nc
   d(i).ps    = dir([d(i).name '*.ps']); % Any Other Plots
   d(i).lad   = [d(i).name '.lad'];
   f = netcdf([loc '/' d(i).nc]);
   dims = dim(f);
   atts = att(f);
   
   % Find All Possible Variables (Sorted By Dimension)
   for j = 1:length(dims)
      field = ncnames(dims{j});
      dum = char(ncnames(var(dims{j}))');
      if ~isempty(dum)  
	 for k = 1:size(dum,1)
	    eval(['H.' dum(k,:) ' = datatype(f{''' deblank(dum(k,:)) '''});']);
	 end
	 if ~isfield(F, field)
	    eval(['F.' char(field) ' = dum;']);
	    eval(['G.' char(field) ' = length(f(''' char(field) '''));']);
	 else
	    eval(['G.' char(field) ' = max(G.' char(field) ', ' ...
		  'length(f(''' char(field) ''')));']);	 	 
	    dum = unique(char([cellstr(dum); ...
			       cellstr(eval(['F.' char(field)]))]), 'rows');
	    eval(['F.' char(field) ' = dum;']);
	 end
      end
   end
   
   % Find Special Global Attributes -> Cast Variables
   for j = 1:length(atts)
      at = char(ncnames(atts{j}));
      if any(strncmp(ncnames(atts{j}), AttVar,4))
	 val = eval(['f.' at '(:);']); 
	 if ischar(val)
	    eval(['F.' at ' = at;'])
	    eval(['H.' at ' = ''char'';']);
	    if ~isfield(G, at)
	       eval(['G.' at ' = length(val);']);
	    else
	       eval(['G.' at ' = max(length(val), G.' at ');']);
	    end
	 elseif numel(val) > 1
	    eval(['F.' at ' = at;']);
	    eval(['H.' at ' = ''float'';'])
	    if ~isfield(G, at)
	       eval(['G.' at ' = length(val);']);
	    else
	       eval(['G.' at ' = max(length(val), G.' at ');']);
	    end
	 else
	    eval(['H.' at ' = ''float'';'])
	    if ~isfield(F, 'cast')
	       eval(['F.cast = at;' ]);
	    else
	       F.cast = unique(char([cellstr(at); cellstr(F.cast)]), 'rows');
	    end
	 end
      end
   end
   close(f);
end


% Add Interpolation Grids 
field = fieldnames(GRID);
for i = 1:length(field)
   eval(['G.' field{i} ' = length(GRID.' field{i} ');']);
end

% Merge Dimensions of Length 1 -> 'cast'
field = fieldnames(G);
dd = [];
for i = 1:length(field)
   if eval(['G.' field{i} ' == 1']) % Dimension Length = 1
      dd = [dd i];
   end
end
if ~isempty(dd)
   ff = field(dd);
   if isfield(F, 'cast')
      dd = cellstr(F.cast);
   else
      dd = [];
   end
   for i = 1:length(ff)
      dd = [dd; cellstr(eval(['F.' ff{i}]))];
   end
   for i = 1:length(ff)
      eval(['F = rmfield(F, ''' ff{i} ''');']);
      eval(['G = rmfield(G, ''' ff{i} ''');']);
   end
   F.cast = char(dd);
   G.cast = 1;
end

dd = fieldnames(G);
disp([' Found ' int2str(length(dd)) ' Dimensions'])
%for i = 1:length(dd)
%   disp(['  ' dd{i}]);
%end

dd = fieldnames(H);
disp([' Found ' int2str(length(dd)) ' Variables'])
%for i = 1:length(dd)
%   disp(['  ' dd{i}]);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Dummy Arrays For Every Cast %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' Creating Dummy Arrays For Every Cast');
d1 = length(d);
dims = fieldnames(F);
for i = 1:length(dims)
   if isfield(GRID, dims{i})
      d2 = length(eval(['GRID.' dims{i} ]));
   else
      d2 = eval(['G.' dims{i}]);
   end
   vars = eval(['F.' dims{i}]);
   if ~isempty(vars)
      for j = 1:size(vars,1)
	 if eval(['strcmp(H.' vars(j,:) ', ''char'')'])
	    eval(['D.' vars(j,:) '= repmat('' '', d1, d2);']);	    
	 else
	    eval(['D.' vars(j,:) '= repmat(nan, d1, d2);']);
	 end
      end
   else
      eval(['F = rmfield(F, ''' vars{i} ''');']);
   end
end
field = fieldnames(GRID);
for i = 1:length(field)
   eval(['D.' field{i} '=GRID.' field{i} ';']);
end
D.cast = [1:length(d)]';

%%%%%%%%%%%%%%%%%%
% Load Cast Data %
%%%%%%%%%%%%%%%%%%

disp(' Loading Cast Data');
for i = 1:length(d)
   fprintf(['  ' d(i).nc]);
   f = netcdf([loc '/' d(i).nc]);
   dims = ncnames(dim(f));
   vars = ncnames(var(f));
   atts = ncnames(att(f));
   D.lat(i,1) = f.lat(:);
   D.lon(i,1) = f.lon(:);
   D.maxz(i,1) = max(f{'z'}(:));

   % Global Attributes -> Variables
   for j = 1:length(atts)
      if any(strncmp(atts{j}, {'GEN_'; 'BAR_'; 'INPU'; 'LADC'},4))
	 nam = atts{j};
	 val = eval(['f.' nam '(:)']);
	 eval(['D.' nam '(i,1:length(val)) = val;']);
      end      
   end 

   % Catalog Variables
   D.cast_id(i,1:length(f{'name'}(:))) = f{'name'}(:)';
   D.station(i,1) = f.ladcp_station(:);
   start = julian(2000,1,1,0);
   D.date(i,:) = f{'date'}(:)';
   D.jul(i,1)  = julian(D.date(i,:)) - start;

   % Additional Variables
   for j = 1:length(dims)
     dd = dims{j};
     if isfield(F, dd)
       v = cellstr(eval(['F.' dd]));
       for k = 1:length(v)
	 vv = v{k};
	 %clearline
	 %fprintf(['  ' d(i).nc ' : ' dd ' : ' vv ])
	 if any(strcmp(vv, vars)) & ~strcmp(vv,dd) % Found a Variable
	   if isfield(GRID, dd)  % We need to Interpolate
	     oldgrid = eval(['f{''' dd '''}(:)']);
	     switch dd
	      case 'zbot'
	       oldgrid = eval(['f.GEN_Ocean_depth_m(:)']) - oldgrid;
	      case 'tim'
	       oldgrid = (oldgrid - oldgrid(1)) / (oldgrid(end)-oldgrid(1));
	     end
	     olddata = eval(['f{''' vv '''}(:)']);
	     newgrid = eval(['GRID.' dd]);
	     eval(['D.' vv '(i,:) = interp1q(oldgrid, olddata, newgrid);']);
	   else % Just Store the Data
	     data = eval(['f{''' vv '''}(:)']);
	     ii = length(data);
	     eval(['D.' vv '(i,1:ii) = data;']);
	   end
	 end
	 if any(strcmp(vv, vars)) % Variable or Dimension
           % Store Attributes
	   aa = eval(['ncnames(att(f{''' vv '''}));']);   
	   for t = 1:length(aa)
	     a = aa{t};
	     eval(['A.' vv '.' a ' = f{''' vv '''}.' a '(:);']);
	   end
	   % Non-char Variables : Add missing_value Tag
	   if eval(['~strcmp(H.' vv ', ''char'') & ~strcmp(vv,dd)']) 
	      eval(['A.' vv '.missing_value = nan;']);
	   end
	 end
       end
%     end
     %%%  
     elseif isfield(H,dd)
        % Store Attributes
	aa = eval(['ncnames(att(f{''' dd '''}));']);   
	for t = 1:length(aa)
	   a = aa{t};
	   eval(['A.' dd '.' a ' = f{''' dd '''}.' a '(:);']);
	end
        % Non-char Variables : Add missing_value Tag
	if eval(['~strcmp(H.' dd ', ''char'')']) 
	   eval(['A.' dd '.missing_value = nan;']);
	end
     end
     %%%
   end
   close(f)
   clearline
   fprintf(['  ' d(i).nc '\n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Attributes For Variables Created Here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = rmfield(F, 'name'); %
G = rmfield(G, 'name'); % For Ingrid Compatability
H = rmfield(H, 'name'); % We have renamed 'name' to 'cast_id'
A = rmfield(A, 'name'); %

F.cast     = char([cellstr(F.cast); cellstr({'jul';'maxz';'station'})]);
F.cast_id  = 'cast_id';
G.cast     = size(D.cast,1);
G.cast_id  = size(D.cast_id,2);
H.cast     = 'float';
H.cast_id  = 'char';
H.jul      = 'float';
H.maxz     = 'float';
H.station  = 'float';
A.cast.long_name    = 'Cast';
A.cast_id.long_name = 'Cast ID';
A.jul.long_name     = 'Date';
A.jul.units         = 'Days Since 2000-01-01';
A.maxz.long_name    = 'Maximum Depth';	
A.maxz.units        = 'Meters';
A.station.long_name = 'Station Number';

% Limits for Ingrid Plot
fac = 0.1;
lon = D.lon;
lon(lon < 0) = lon(lon < 0) + 360;
lonmin = floor(min(lon));
lonmax = ceil(max(lon));
lonmin = lonmin - fac*(lonmax-lonmin);
lonmax = lonmax + fac*(lonmax-lonmin);
latmin = floor(min(D.lat));
latmax = ceil(max(D.lat));
latmin = latmin - fac*(latmax-latmin);
latmax = latmax + fac*(latmax-latmin);
A.lon.scale_min = lonmin;
A.lon.scale_max = lonmax;
A.lat.scale_min = latmin;
A.lat.scale_max = latmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Cruise netCDF File %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' Creating Cruise NetCDF matlab');
eval([' save ',loc '/' name ' A D F G H ']);

disp(' Creating Cruise NetCDF File');
f = netcdf([loc '/' name '.nc'], 'clobber');

% Dimensions
disp(' Saving Dimensions');
FF = fieldnames(F);
for i = 1:length(FF) 
   dd = FF{i};
   DD = dd;
   % Convert BLAH_blah to BLAH.blah if BLAH is in AttVar
   % This is done to create subgroups for Ingrid
   if any(strncmp(DD, AttVar,4))
      ii = min(find(DD == '_'));
      if ~isempty(ii)
	 DD(ii) = '.';
      end
   end
   SKIP = 0;
   if eval(['strcmp(H.' dd ', ''char'')'])
      DD = [DD '_char'];
      SKIP = 1;
   elseif eval(['numel(D.' dd ') ~= length(D.' dd ')'])
      DD = [DD '_num'];      
      SKIP = 1;
   end
   eval(['f(''' DD ''') = G.' dd ';']);
   if ~SKIP
      eval(['f{''' DD '''} = ncfloat(''' DD ''');']);
   end
   if eval(['numel(D.' dd ') == length(D.' dd ')'])
      eval(['f{''' DD '''}(:) = D.' dd ';']);      
   else
% Fixes ingrid problem      
%      eval(['f{''' DD '''}(:) = [1:G.' dd '];']);
   end
   % Dimension Attributes
   if isfield(A, dd)
      aa = eval(['fieldnames(A.' dd ')']);
      for j = 1:length(aa)
	 at = aa{j};
	 eval(['f{''' DD '''}.' at '=A.' dd '.' at ';']);
      end
   end
   %disp(['  ' DD]);
end

% Variables
disp(' Saving Variables');
FF = fieldnames(F);
for i = 1:length(FF) 
   dd = FF{i};
   DD = dd;
   if eval(['strcmp(H.' dd ', ''char'')'])
      DD = [DD '_char'];
   elseif eval(['numel(D.' dd ') ~= length(D.' dd ')'])
      DD = [DD '_num'];      
   end
   if any(strncmp(DD, AttVar,4))
      ii = min(find(DD == '_'));
      if ~isempty(ii)
	 DD(ii) = '.';
      end
   end
   V = eval(['F.' dd]);
   for j = 1:size(V,1)
      vv = char(cellstr(V(j,:)));
      VV = vv;
      % What Kind of Variable?
      % Main: z,u,v, etc.
      % Extra: u_do, range, nvel, etc.
      % Attribute Derived: GEN_Depth_source, BAR_ref_descr, etc.

      if ~strcmp(DD, vv) % Variable, Not a Dimension
         % Attribute Derived
	 if any(strncmp(VV, AttVar,4))
	    ii = min(find(VV == '_'));
	    if ~isempty(ii)
	       VV(ii) = '.';
	    end
         % Extra
	 elseif ~any(strcmp(VV, MainVariables))
	    switch dd
	     case 'z'
	      VV = ['EXTRA.PROFILE.' VV];
	     case 'tim'
	      VV = ['EXTRA.TIMESERIES.' VV];
	     case 'zbot'
	      VV = ['EXTRA.BOTTOM.' VV];
	     case 'z_sadcp'
	      VV = ['EXTRA.SADCP.' VV];
	    end
	 end 
	 if ~strcmp(vv, 'cast') & ~strcmp(dd, 'cast')
	    if eval(['strcmp(H.' vv ', ''char'')'])
	       eval(['f{''' VV '''} = ncchar(''cast'', ''' DD ''');']);
	    else
	       eval(['f{''' VV '''} = ncfloat(''cast'', ''' DD ''');']);
	    end
	    eval(['f{''' VV '''}(:,:) = D.' vv ';']);
	 else
	    if eval(['strcmp(H.' vv ', ''char'')'])
	       eval(['f{''' VV '''} = ncchar(''cast'');']);
	    else
	       eval(['f{''' VV '''} = ncfloat(''cast'');']);	    
	    end
	    eval(['f{''' VV '''}(:) = D.' vv ';']);
	 end
         % Variable Attributes
	 if isfield(A, vv)
	    aa = eval(['fieldnames(A.' vv ')']);
	    for k = 1:length(aa)
	       at = aa{k};
	       eval(['f{''' VV '''}.' at '=A.' vv '.' at ';']);
	    end
	 end
      else

      end
      %disp(['  ' VV]);
   end
end

% Global Attributes
f.EXTRA.description            = 'More Data';
f.EXTRA.PROFILE.description    = 'PROFILES';
f.EXTRA.TIMESERIES.description = 'TIMESERIES';
f.EXTRA.BOTTOM.description     = 'BOTTOM';
f.BAR.description              = 'Barotropic Data';
f.GEN.description              = 'General Processing Parameters';
f.INPUT.description            = 'Available Inputs';
f.LADCP.description            = 'LADCP Technical Data';

close(f);

% Create TAR Files for .nc, .lad, .ps files
disp(' Finished')
htmltable(directory, dods)

function [] = clearline()
fprintf(repmat('\b',1,80));
fprintf(repmat(' ' ,1,80));
fprintf(repmat('\b',1,80));

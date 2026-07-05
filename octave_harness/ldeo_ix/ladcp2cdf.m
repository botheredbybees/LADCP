function [] = ladcp2cdf(fname,dr_struct,da,p,ps,f,att);
% function [] = ladcp2cdf(fname,dr_struct,da,p,ps,f,att);
%
% function to save LADCP data into a netcdf file for MatLab version 2012a
%
% input  :	fname		- output filename
%			dr_struct	- main inversion results (velocity profiles)
%			da,p,ps,f,att - arbitrary metadata structures
%
% Subroutine :  add_struct

% Created By:   Diana Cardoso, Bedford Institute of Oceangraphy
%               Diana.Cardoso@dfo-mpo.gc.ca
% Description:  Based on LDEO software to Process LADCP, version IX.8,
%               script ladcp2cdf.m version 0.1	last change 08.03.2002. 
%               maintained by A.M. Thurnherr and downloaded from:
%       http://www.ldeo.columbia.edu/cgi-bin/ladcp-cgi-bin/hgwebdir.cgi
%       The function ladcp2cdf was changed to run with the the Matlab
%       version 2012, which now supports netcdf.

%======================================================================
%                    L A D C P 2 C D F . M 
%                    doc: Thu Aug 15 10:52:55 2013
%                    dlm: Thu Nov 26 22:01:32 2015
%                    (c) 2013 A.M. Thurnherr, from code contributed by D. Cardoso
%                    uE-Info: 183 0 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% NOTES:
%	- This version creates slightly different files than the original version
%	  created by Visbeck/Krahmann. In the original version, the contents of the
%	  dr structure end up as top-level variables and the contents of
%	  the da, p, ps, f and att structures end of as global attributes. In 
%	  the new version, the latter are saved as sub-structures, with _struct appended
%	  to the internal names to avoid conflicts.

% HISTORY:
%   Aug 15, 2013: - incorporated this code, supplied Diana Cardoso, into IX_10beta
%		  - modified doc in header
%		  - renamded struct variable to dr_struct
%		  - removed 'cd' in and out of results directory (pathnames work just fine)
%		  - delete netcdfile before it is written to (old 'clobber' option)
%		  - removed 'l' suffix from all dims
%		  - replaced yes/no logical vals by true/false
%		  - renamed substructures from st2..st6 to internal names (da,p,ps,f,att)
%   Aug 28, 2013: - incorporated bug fix provided by Diana Cardoso to prevent lat,lon,name and 
%		    date to be stored 2cd in the nc file, which can make the code
%		    bomb if the length of any other var is 6 (or equal to the length of name?)
%   Nov     2015: - new coded provided by Diana Cardoso (with support from Eric Firing)
%   Nov 26, 2015: - BUG: code did not work when bottom-track data were missing
%--------------------------------------------------------------------------
% check arguments and remove existing NetCDF file with fname(output filename)
%--------------------------------------------------------------------------

if nargin<2
  error('need two input arguments')
end
if ~isstruct(dr_struct)
  error('second argument must be a dr structure')
end

netcdfile = deblank(fname); %remove any blanks from string end
if exist(netcdfile,'file')
	delete(netcdfile)
end
%--------------------------------------------------------------------------
%Create a classic format NetCDF file with 8 dimension definitions.
%--------------------------------------------------------------------------

mySchema.Name   = '/'; % indicating the full file as opposed to a group
mySchema.Format = 'classic'; %The format of the NetCDF file
%Create Dimensions
mySchema.Dimensions(1).Name   = 'name';
mySchema.Dimensions(1).Length = length(getfield(dr_struct,'name'));
mySchema.Dimensions(2).Name   = 'date';
mySchema.Dimensions(2).Length = 6;
mySchema.Dimensions(3).Name   = 'lat';
mySchema.Dimensions(3).Length = 1;
mySchema.Dimensions(4).Name   = 'lon';
mySchema.Dimensions(4).Length = 1;

if isfield(dr_struct,'z');
  lz = length(getfield(dr_struct,'z'));
else
  lz = 0;
end  
mySchema.Dimensions(5).Name   = 'z';
mySchema.Dimensions(5).Length = lz;

if isfield(dr_struct,'tim');
  ltim = length(getfield(dr_struct,'tim'));
else
  ltim = 0;
end  
mySchema.Dimensions(6).Name   = 'tim';
mySchema.Dimensions(6).Length = ltim;

if isfield(dr_struct,'zbot');
  lbot = length(getfield(dr_struct,'zbot'));
else
  lbot = 0;
end
mySchema.Dimensions(7).Name   = 'zbot';
mySchema.Dimensions(7).Length = lbot;

if isfield(dr_struct,'z_sadcp');
  lsadcp = length(getfield(dr_struct,'z_sadcp'));
else
  lsadcp = 0;
end
mySchema.Dimensions(8).Name   = 'z_sadcp';
mySchema.Dimensions(8).Length = lsadcp;

%--------------------------------------------------------------------------
%Add first 4 variables definitions(name, date, lat, lon) to Variable field
%--------------------------------------------------------------------------

mySchema.Variables(1).Name = 'name';
mySchema.Variables(1).Dimensions(1) = mySchema.Dimensions(1);
mySchema.Variables(1).Datatype = 'char';
mySchema.Variables(2).Name = 'date';
mySchema.Variables(2).Dimensions(1) = mySchema.Dimensions(2);
mySchema.Variables(2).Datatype = 'int32';
mySchema.Variables(3).Name = 'lat';
mySchema.Variables(3).Dimensions(1) = mySchema.Dimensions(3);
mySchema.Variables(3).Datatype = 'single';
mySchema.Variables(4).Name = 'lon';
mySchema.Variables(4).Dimensions(1) = mySchema.Dimensions(4);
mySchema.Variables(4).Datatype = 'single';

dr_struct = orderfields(dr_struct) ; %sort dr_struct by fieldnames

fnames = fieldnames(dr_struct);
% find name, date, lat, lon in fnames that have already been added to the 
% NetCDF file and remove from fnames
nn=strncmp('name',fnames,6);   
nda=strncmp('date',fnames,4);
nla=strncmp('lat',fnames,3); 	
nlo=strncmp('lon',fnames,3);
ntot=[nn+nda+nla+nlo]; Ktot = logical(ntot);
fnames(Ktot,:)=[];	

%--------------------------------------------------------------------------
%Add the remaining variables definitions from dr_struct to Variable field
%--------------------------------------------------------------------------

for n=1:size(fnames,1)
    dummy = getfield(dr_struct,fnames{n});
    mySchema.Variables(4+n).Name = fnames{n};
    if ~isempty(strfind(fnames{n},'bot'))
        mySchema.Variables(4+n).Dimensions(1) = mySchema.Dimensions(7);
    elseif ~isempty(strfind(fnames{n},'sadcp'))
        mySchema.Variables(4+n).Dimensions(1) = mySchema.Dimensions(8);
    elseif ~isempty(strfind(fnames{n},'tim')) ||  ~isempty(strfind(fnames{n},'ship')) ...
        || ~isempty(strfind(fnames{n},'ctd')) && isempty(strmatch('ctd_t',fnames{n},'exact')) ...
        && isempty(strmatch('ctd_s',fnames{n}))
        mySchema.Variables(4+n).Dimensions(1) = mySchema.Dimensions(6);
    elseif ~isempty(strfind(fnames{n},'vel')) ||  ~isempty(strfind(fnames{n},'shear')) ...
        || ~isempty(strfind(fnames{n},'ctd_')) || ~isempty(strfind(fnames{n},'uerr')) ...
        ||  ~isempty(strfind(fnames{n},'range'))  || ~isempty(strfind(fnames{n},'ts')) ...
        || ~isempty(strmatch('p',fnames{n},'exact')) || ~isempty(strmatch('u',fnames{n},'exact')) ...
        || ~isempty(strmatch('v',fnames{n},'exact'))|| ~isempty(strmatch('z',fnames{n},'exact')) ...
        || ~isempty(strfind(fnames{n},'u_')) || ~isempty(strfind(fnames{n},'v_'));
        mySchema.Variables(4+n).Dimensions(1) = mySchema.Dimensions(5);
    else
        mySchema.Variables(4+n).Dimensions.Name = fnames{n};
        mySchema.Variables(4+n).Dimensions.Length = length(dummy);
    end 
    if strncmp(fnames{n},'tim',4)==1 ||  ~isempty(strfind(fnames{n},'vbar')) ||  ~isempty(strfind(fnames{n},'ubar'))
        mySchema.Variables(4+n).Datatype = 'double';
    elseif strncmp(fnames{n},'nvel',4)==1 
        mySchema.Variables(4+n).Datatype = 'int16';
    else
        mySchema.Variables(4+n).Datatype = 'single'; 
    end   
%    disp(sprintf('%d: %s %d',n,mySchema.Variables(4+n).Dimensions.Name,mySchema.Variables(4+n).Dimensions.Length));
end

%--------------------------------------------------------------------------
%Add the attributes definitions from dr_struct to the Variable field
%--------------------------------------------------------------------------

ncwriteschema(fname, mySchema);
ncwriteatt(fname,'name','long_name',att.name.long_name);
ncwriteatt(fname,'date','long_name',att.date.long_name);
ncwriteatt(fname,'date','units',att.date.units);
ncwriteatt(fname,'lat','long_name',att.lat.long_name);
ncwriteatt(fname,'lat','units',att.lat.units);
ncwriteatt(fname,'lon','long_name',att.lon.long_name);
ncwriteatt(fname,'lon','units',att.lon.units);
ncwrite(fname,'name',dr_struct.name);
ncwrite(fname,'date',dr_struct.date);
ncwrite(fname,'lat',dr_struct.lat);
ncwrite(fname,'lon',dr_struct.lon);

for n=1:size(fnames,1)
  dummy = getfield(dr_struct,fnames{n});
    ncwrite(fname,fnames{n},dummy);
    if isfield(att,fnames{n});
        ncwriteatt(fname,fnames{n},'long_name',att.(fnames{n}).long_name);
        if isfield(att.(fnames{n}),'units');
             ncwriteatt(fname,fnames{n},'units',att.(fnames{n}).units);
        end
    end
end

%--------------------------------------------------------------------------
%Add the attributes definitions from da, p, f, ps mat structures to
%Attributes field using subroutine add_struct
%--------------------------------------------------------------------------

% remove duplicate fields from mat structures
ps = rmfield(ps,'outlier');
p = rmfield(p,'checkpoints');
f = rmfield(f,'ctd_time_base');
f = rmfield(f,'nav_time_base');

% combine all mat structures and sort fieldnames
names = [fieldnames(da); fieldnames(p)];
struct1 = cell2struct([struct2cell(da); struct2cell(p)], names, 1);
names = [fieldnames(struct1); fieldnames(ps)];
struct2 = cell2struct([struct2cell(struct1); struct2cell(ps)], names, 1);
names = [fieldnames(struct2); fieldnames(f)];
struct3 = cell2struct([struct2cell(struct2); struct2cell(f)], names, 1);
structsorted = orderfields(struct3) ;

% add attributes to netcdf from the p, ps, da and f structures
% the slash indicates a global variable, if you change it the value in the 
% structures will be placed in the Variables field of the netcdf 
add_struct(fname,'/',structsorted); 

end % function

%----------------------------------------------------------------------

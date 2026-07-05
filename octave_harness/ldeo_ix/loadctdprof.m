%======================================================================
%                    L O A D C T D P R O F . M 
%                    doc: Thu Jun 17 23:17:25 2004
%                    dlm: Fri Mar  5 15:48:54 2010
%                    (c) 2004 ladcp@
%                    uE-Info: 77 13 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jul 17, 2008: - removed dependence on p.nav_start, because that's
%		    not set in [loadnav.m] any more
%   Jul 27, 2008: - nanmean() -> meannan()
%   Jan  7, 2009: - tightened use of exist()

function [d,p]=loadctdprof(f,d,p,ipos,isub)
% function [d,p]=loadctdprof(f,d,p,ipos,isub)
% LADCP-2 software version 7.0 
%
% ========== THIS PART IS SPECIFIC TO LDEO CTD PROFILE DATA ============
%
% you might want to change it to accomodate your own data format
%
%
% load and merges LDEO ctdfile
% ipos(1) : pressure
% ipos(2) : temperature (in situ)
% ipos(3) : salinity
%
% computes sound speed 
% computes static stability
%
% Martin Visbeck, 6/10/99
% revised March 2002
% revised December 2002

% CLIVAR P02 (VANC32) CTD profiles do not contain lat/lon;
% threfore, they are taken from p.nav_start

if nargin<5, isub=1; end
if nargin<4, ipos=[1,2,3]; end

% read ctd data file
disp(['LOADCTDPROF: load CTD profile ',f.ctdprof])
if ~exist(f.ctdprof,'file')
 warn=([' LOADCTDPROF can not find ',f.ctdprof]);
 p.warnp(size(p.warnp,1)+1,1:length(warn))=warn;
 disp(warn)
 return
end
[dctd,posctd]=readctd(f.ctdprof,isub,ipos);

if length(posctd) == 2
	d.ctdprof_lon=posctd(2);
	d.ctdprof_lat=posctd(1);
elseif p.navdata
	d.ctdprof_lon = meannan(d.slon);
	d.ctdprof_lat = meannan(d.slat);
else
	error('do not know how to determine CTD position');
end
d.ctdprof_z=p2z(dctd(:,1),d.ctdprof_lat);
d.ctdprof_p=dctd(:,1);
d.ctdprof_t=dctd(:,2);
d.ctdprof_s=dctd(:,3);

% get N^2
if exist('sw_bfrq','file')
 d.ctdprof_N2=sw_bfrq(d.ctdprof_s,d.ctdprof_t,d.ctdprof_p,d.ctdprof_lat);
 d.ctdprof_N2(end+1)=d.ctdprof_N2(end);
 d.ctdprof_ss=sw_svel(d.ctdprof_s,d.ctdprof_t,d.ctdprof_p);
else
 disp(' download SW routines to get soundspeed and N^2 ')
end

% get maximum ctdepth
disp([' CTD max depth : ',int2str(max(d.ctdprof_z))])
if ~isfinite(p.zpar(2));
 p.zpar(2)=fix(max(d.ctdprof_z));
end


%====================================================================

function depth=p2z(p,lat)
% !!!!!! USES Z=0 AT P=0  (I.E. NOT 1ATM AT SEA SURFACE)
%	pressure to depth conversion using
%	saunders&fofonoff's method (deep sea res.,
%	1976,23,109-111)
%	formula refitted for alpha(p,t,s) = eos80
%	units:
%		depth         z        meter
%		pressure      p        dbars  (original in bars, but below
%                                              division by 10 is included)
%		latitude      lat      deg
%	checkvalue:
%		depth =       9712.654  m
%	for
%		p     =         1000.     bars
%		lat   =           30.     deg
%	real lat,p
        if nargin < 2, lat=54; end
        p=p/10.;
	x=sin(lat/57.29578);
	x=x*x;
	gr=9.780318*(1.0+(5.2788e-3+2.36e-5*x)*x)+1.092e-5*p;
	depth=(((-1.82e-11*p+2.279e-7).*p-2.2512e-3).*p+97.2659).*p;
	depth=depth./gr;
%
%
%====================================================================

function [d,pos]=readctd(file,isub,index)
% read CTD data
id=fopen(file,'r');
A=fread(id);
fclose(id);
A=setstr2(A);

pos = [];

%%%ii=find(A(:,1)=='@');
%%%A=A((ii+1):end,:);

[lt,lv]=size(A);

d=[];

if nargin<2, isub=1; end

for i=1:isub:lt
 d=[d;sscanf(A(i,:),'%g')'];
end

if nargin>2
 d=d(:,index);
end


% -----------------------------------------------------
function        AR=setstr2(A,cr)
% function      AR=setstr2(A,cr)
% reshape charater string vector A
% in record structure
% cr = 10 (default)

% remove line feeds
ii=find(A==13);
A(ii)=[];


if nargin<2, cr=10;, end

ii=[0;find(A==cr)];

n=length(ii)-1;
m=max([diff(ii)]);

% fill array with blanks
AR=ones(n,m)*32;

for i=1:n
 j=(ii(i)+1):(ii(i+1)-1);
 if length(j)>0
  j1=j-j(1)+1;
  AR(i,j1)=A(j)';
 end
end

AR=setstr(AR);


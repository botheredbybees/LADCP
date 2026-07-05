%======================================================================
%                    H T M L T A B L E . M 
%                    doc: Wed Jun  4 12:58:37 2008
%                    dlm: Wed Nov 18 20:44:02 2009
%                    (c) 2008 M. Visbeck & A.M. Thurnherr
%                    uE-Info: 49 73 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

function [] = htmltable(directory, dods)
% function [] = htmltable(directory, dods)
%
% Create HTML Table of LADCP Files in a Directory
%
% Input :  directory  : Either a directory '/home/someone/lookhere/'
%                       Or the 'f' structure from LADCP where
%                       'f.res' will be used.
%      
%          dods       : Structure with DODS info ** OPTIONAL **
%                       dods.cruise - dods address for cruise file
%                       dods.d1     - start of address for cast
%                       dods.d2     - end of address for cast
%                       
% Output : Files created in your directory:
%             index.html  - Start Page with Map
%             index.htm   - Start Page with Map
%             table0.html - Start Page with Map
%             table1.html - Individual Cast Data
%             table2.html - Processing Methods
%             table3.html - Technical
%             table4.html - Plots
%             cruise.ps   - Map of Cruise (filename will vary)
%             cruise.gif  - Map of Cruise (filename will vary)
%
% For DODS addresses, htmltable will insert the cast number to get the
% correct address.
%
% Example: directory = '../ladcp_output/';
%
%          dods.cruise = ['http://kage.ldeo.columbia.edu:81/SOURCES/' ...
%                         '.LDEO/.ClimateGroup/.LADCP/.test/dods'];
%          dods.d1     = ['http://kage.ldeo.columbia.edu:81/SOURCES/' ...
%                         '.LDEO/.ClimateGroup/.LADCP/.test/N+'];
%          dods.d2     = ['+VALUE/dods'];
%

% Modifications by ANT:
%	Jun  4, 2008: - replace all spaces by _ in cruise name
%   Jan  7, 2009: - tightened use of exist()
%	Nov 18, 2009: - replaced "convert" by "pstogif" (syntax is identical)

if nargin < 1, help htmltable, return, end
if nargin == 2 & isstruct(dods)
   DODS = 1;
else
   DODS = 0;
end

disp('HTMLTABLE')
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

% Ignore Cruise File
count = 0;
bad = [];
for i = 1:length(d)
   f = netcdf([loc '/' d(i).name]);
   % Is this an LADCP File?
   if isempty(f.GEN_Profile_start_decimal_day(:))
      bad = [bad i];
      if ~isempty(f{'cast'}(:))
	 cruise_cast = f{'cast'}(:);
	 cruise_station = f{'station'}(:);
      end
   else
      count = count + 1;
   end
   close(f)
end
d(bad) = [];
if ~exist('cruise_cast','var')
   DODS = 0;
end

if count 
   disp([' Found ' int2str(count) ' Files to Process']);
   for i = 1:length(d)
      disp(['  ' d(i).name]);
   end
else
   error([' Nothing To Do!';])
end
f = netcdf([loc '/' d(1).name]);
name = regexprep(f.cruise_id(:),' ','_');
soft = f.software(:);
time = date;
time(time == '-') = ' ';
if isempty(name)
   name = 'unknown';
end
disp([' Cruise ID: ' name ]);
disp([' Software:  ' soft ]);
disp([' Creating HTML Tables']);

% Header
cr=char(10); % Carriage Return

H0 = ['<TITLE>' name '</TITLE>' cr ...
      '<B>' name '</B>' cr '<BR>' cr ...
      'Created on ' time ' using ' soft cr '<BR>' ...
      '<a href="table0.html"> Cruise Data </a> &nbsp' cr ...
      '<a href="table1.html"> Individual Cast Data </a> &nbsp' cr ...
      '<a href="table2.html"> Processing Methods</a> &nbsp' cr ...
      '<a href="table3.html"> Technical </a> &nbsp' cr ...
      '<a href="table4.html"> Plots </a> &nbsp' cr ...
      '<a href="../index.html"> Other Cruises </a> &nbsp' cr];

% TABLE 0
H1 = [H0, '<img SRC="' name '.gif" HEIGHT=800 WIDTH=600 ALIGN=center>' cr];
if DODS
   H1 = [H1 '<br><ul><li>' cr...
	 '<a href="' dods.cruise(1:end-4) ...
	 '">Live Access IRI/LDEO Data Server</a><br>'];
else
      H1 = [H1 '<br><ul><li>' cr...
	 'Live Access IRI/LDEO Data Server<br>'];
end

H1 = [H1 '<li>All Cruise Data in one file: ' cr];
H1 = [H1 '<a href="' name '.nc"> NetCDF </a>'];
H1 = [H1 ' <br>' cr '</li><li>' cr];
if DODS
   H1 = [H1 ' <a href="' dods.cruise '"> DODS </a> <br>' cr];   
else
   H1 = [H1 ' DODS <br>' cr];
end
H1 = [H1 '</li><li>' ...
      'One File Per Cast: '];
H1 = [H1 '<a href="' name '_nc.tgz"> NetCDF tarball </a>'];
H1 = [H1 '</li><li>' ];
H1 = [H1 '<a href="' name '_ps.tgz"> PS tarball </a>'];
H1 = [H1 '</li><li>' ];
H1 = [H1 '<a href="' name '_log.tgz"> Logfile tarball </a>'];
H1 = [H1 '</li><li>' ];
H1 = [H1 '<a href="' name '_lad.tgz"> ASCII results tarball </a>'];
H1 = [H1 '</li><li>' ];
H1 = [H1 '<a href="' name '_mat.tgz"> MAT results tarball </a>'];
H1 = [H1 '</li><li>' ];
H1 = [H1 'To unpack tarballs :   tar -xzvf filename.tgz'];

cwd = pwd;
cd(loc)
disp(['  NETCDF tarball']);
system(['tar -czf ' name '_nc.tgz *.nc']); 
disp(['  PostScript tarball']);
system(['tar -czf ' name '_ps.tgz *.ps']); 
disp(['  Logfile tarball']);
system(['tar -czf ' name '_log.tgz *.log']); 
disp(['  ASCII tarball']);
system(['tar -czf ' name '_lad.tgz *.lad']); 
disp(['  Mat-file tarball']);
system(['tar -czf ' name '_mat.tgz *.mat']); 
cd(cwd)

fid = fopen([loc '/table0.html'], 'w');
fwrite(fid, H1);
fclose(fid);
fid = fopen([loc '/index.html'], 'w');
fwrite(fid, H1);
fclose(fid);
fid = fopen([loc '/index.htm'], 'w');
fwrite(fid, H1);
fclose(fid);
disp(['  index.html']);
disp(['  index.htm']);
disp(['  table0.html']);

% TABLE 1
H1 = [H0 '&nbsp ',cr,...
      '<div ALIGN=right><table BORDER COLS=12 WIDTH="100%" NOSAVE >',cr];
H1 = [H1 '<tr align=center>' cr];
H1 = [H1 '<td> Name           </td>' cr];
H1 = [H1 '<td> Year           </td>' cr];
H1 = [H1 '<td> Month          </td>' cr];
H1 = [H1 '<td> Day            </td>' cr];
H1 = [H1 '<td> Lat            </td>' cr];
H1 = [H1 '<td> Lon            </td>' cr];
H1 = [H1 '<td> Max Z          </td>' cr];
H1 = [H1 '<td> Median U Error </td>' cr];
H1 = [H1 '<td> NETCDF File    </td>' cr];
H1 = [H1 '<td> Matlab File    </td>' cr];
H1 = [H1 '<td> Text File      </td>' cr];
H1 = [H1 '<td> DODS           </td>' cr];

for i = 1:length(d)
   f = netcdf([loc '/' d(i).name]);
   warn = f.warnings(:);
   WARN(i) = 0;
   % Are there warnings? If so, we will highlight these
   if length(find(warn == char(10))) > 1
      WARN(i) = 1;
   end
   cruiseid = f{'name'}(:)';
   gtime = f{'date'}(:);
   uerr = f{'uerr'}(:);
   lat = f{'lat'}(:);
   lon = f{'lon'}(:);
   station(i) = f.ladcp_station(:);
   LAT(i) = lat;
   LON(i) = lon;
   maxz = max(f{'z'}(:));
   med_err_vel = median(uerr(~isnan(uerr(:))));
   H1 = [H1 '<tr align=center>' cr];
   if WARN(i)
      H1 = [H1 '<td> <font color="red">' cruiseid '</font></td>' cr];      
   else
      H1 = [H1 '<td>' cruiseid                     '</td>' cr];
   end
   H1 = [H1 '<td>' num2str(gtime(1))            '</td>' cr];
   H1 = [H1 '<td>' num2str(gtime(2))            '</td>' cr];
   H1 = [H1 '<td>' num2str(gtime(3))            '</td>' cr];   
   H1 = [H1 '<td>' num3str(lat, 0,2,'0')        '</td>' cr];
   H1 = [H1 '<td>' num3str(lon, 0,2,'0')        '</td>' cr];
   H1 = [H1 '<td>' num2str(maxz)                '</td>' cr];
   H1 = [H1 '<td>' num3str(med_err_vel,0,3,'0') '</td>' cr];   
   H1 = [H1 '<td> <a href="' d(i).name    '"> nc </td>' cr];
   mat = [d(i).name(1:end-2) 'mat'];
   if exist([loc '/' mat],'file')
      H1 = [H1 '<td> <a href="' mat '"> mat </td>' cr];
   else
      H1 = [H1 '<td> mat </td>' cr];
   end
   lad = [d(i).name(1:end-2) 'lad'];
   if exist([loc '/' lad],'file')
      H1 = [H1 '<td> <a href="' lad '"> lad </td>' cr];
   else
      H1 = [H1 '<td> lad </td>' cr];
   end
   if DODS
      cast = cruise_cast(find(cruise_station == station(i)));
      H1 = [H1 '<td> <a href="' dods.d1 int2str(cast) dods.d2 ...
	    '"> dods </a> <td>' cr];      
   else
      H1 = [H1 '<td> dods <td>' cr];
   end
   close(f)
end

H1 = strrep(H1, '<td></td>', '<td>&nbsp </td>'); % Pad Empty Cells
fid = fopen([loc '/table1.html'], 'w');
fwrite(fid, H1);
fclose(fid);
disp(['  table1.html']);


% TABLE 2 - Processing Data
empty = '<td>&nbsp </td>';
H1 = [H0 '&nbsp ',cr,...
      '<div ALIGN=right><table BORDER COLS=6 WIDTH="100%" NOSAVE >', cr];
H1 = [H1 '<tr align=center>' cr];
H1 = [H1 '<td> Name </td>' cr];
H1 = [H1 '<td> Bar Ref </td>' cr];
H1 = [H1 '<td> Depth Source </td>' cr];
H1 = [H1 '<td> Sound Speed </td>' cr];
H1 = [H1 '<td> Percent 3-Beam </td>' cr];
H1 = [H1 '<td> Log </td>' cr];

for i = 1:length(d)
   H1 = [H1 '<tr align=center>' cr];
   f = netcdf([loc '/' d(i).name]);
   cruiseid = f{'name'}(:)';
   if WARN(i)
      H1 = [H1 '<td> <font color="red">' cruiseid '</font></td>' cr];      
   else
      H1 = [H1 '<td>' cruiseid                     '</td>' cr];
   end
   dum = f.BAR_ref_descr(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' dum '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.GEN_Depth_source(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' dum '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.GEN_Sound_sp_calc(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' dum '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.GEN_Percent_3beam(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num2str(dum) '</td>' cr];
   else
      H1 = [H1 empty];
   end
   H1 = [H1 '<td> <a href="' d(i).name(1:end-2) 'log"> log </td>' cr];
   close(f)
end
H1 = strrep(H1, '<td></td>', '<td>&nbsp </td>'); % Pad Empty Cells
fid = fopen([loc '/table2.html'], 'w');
fwrite(fid, H1);
fclose(fid);
disp(['  table2.html']);

% TABLE 3 - LADCP Technical Data
H1 = [H0 '&nbsp',cr,...
      '<div ALIGN=right><table BORDER COLS=11 WIDTH="100%" NOSAVE >', cr];
H1 = [H1 '<tr align=center>' cr];
H1 = [H1 '<td> Name </td>' cr];
H1 = [H1 '<td> LADCP Up Hard Type </td>' cr];
H1 = [H1 '<td> LADCP Up Hard SN </td>' cr];
H1 = [H1 '<td> LADCP Up Hard Conf Single Ping Acc </td>' cr];
H1 = [H1 '<td> LADCP Up Conf Bin Len </td>' cr];
H1 = [H1 '<td> LADCP Up Conf Number Pings </td>' cr];
H1 = [H1 '<td> LADCP Dn Hard Type </td>' cr];
H1 = [H1 '<td> LADCP Dn Hard SN </td>' cr];
H1 = [H1 '<td> LADCP Dn Hard Conf Single Ping Acc </td>' cr];
H1 = [H1 '<td> LADCP Dn Conf Bin Len </td>' cr];
H1 = [H1 '<td> LADCP Dn Conf Number Pings </td>' cr];

for i = 1:length(d)
   f = netcdf([loc '/' d(i).name]);
   H1 = [H1 '<tr align=center>' cr];
   cruiseid = f{'name'}(:)';
   if WARN(i)
      H1 = [H1 '<td> <font color="red">' cruiseid '</font></td>' cr];      
   else
      H1 = [H1 '<td>' cruiseid                     '</td>' cr];
   end
   dum = f.LADCP_up_hard_type(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' dum '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_up_hard_SN(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num2str(dum) '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_up_conf_single_ping_acc(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num3str(dum,1,3,'0') '</td>' cr];
   else
      H1 = [H1 empty];
   end
    dum = f.LADCP_up_conf_bin_len_m(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num2str(dum) '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_up_conf_number_pings(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num2str(dum) '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_dn_hard_type(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' dum '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_dn_hard_SN(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num2str(dum) '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_dn_conf_single_ping_acc(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num3str(dum,1,3,'0') '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_dn_conf_bin_len_m(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num2str(dum) '</td>' cr];
   else
      H1 = [H1 empty];
   end
   dum = f.LADCP_dn_conf_number_pings(:);
   if ~isempty(dum)
      H1 = [H1 '<td>' num2str(dum) '</td>' cr];
   else
      H1 = [H1 empty];
   end
   close(f)
end
H1 = strrep(H1, '<td></td>', '<td>&nbsp </td>'); % Pad Empty Cells
fid = fopen([loc '/table3.html'], 'w');
fwrite(fid, H1);
fclose(fid);
disp(['  table3.html']);

% TABLE 4 - Plots
H1 = [H0 '&nbsp ',cr,...
      '<div ALIGN=right><table BORDER COLS=11 WIDTH="100%" NOSAVE >', cr];
H1 = [H1 '<tr align=center>' cr];
H1 = [H1 '<td> Name                     </td>' cr]; %  Name
H1 = [H1 '<td> Summary                  </td>' cr]; %  1
H1 = [H1 '<td> Engineering Data         </td>' cr]; %  2
H1 = [H1 '<td> Data Quality             </td>' cr]; %  3
H1 = [H1 '<td> Depth                    </td>' cr]; %  4
H1 = [H1 '<td> Heading Corrections      </td>' cr]; %  5
H1 = [H1 '<td> Up/Down Differences      </td>' cr]; %  6
H1 = [H1 '<td> CTD Position             </td>' cr]; %  7
H1 = [H1 '<td> Shear                    </td>' cr]; %  8
H1 = [H1 '<td> SADCP U, V               </td>' cr]; %  9
H1 = [H1 '<td> U, V Offsets, Tilt Error </td>' cr]; % 10

for i = 1:length(d)
   H1 = [H1 '<tr align=center>' cr];
   f = netcdf([loc '/' d(i).name]);
   cruiseid = f{'name'}(:)';
   if WARN(i)
      H1 = [H1 '<td> <font color="red">' cruiseid '</font></td>' cr];      
   else
      H1 = [H1 '<td>' cruiseid                     '</td>' cr];
   end
   for j = 1:10
      ps = char([d(i).name(1:end-3) '_' int2str(j) '.ps']);
      if exist([loc '/' ps],'file')
	 H1 = [H1 '<td><a href=' ps '> PS ' int2str(j) '</a></td>' cr];
      else
	 H1 = [H1 '<td> PS ' int2str(j) '</td>' cr];
      end
   end
   close(f)
end
fid = fopen([loc '/table4.html'], 'w');
fwrite(fid, H1);
fclose(fid);
disp(['  table4.html']);

% Map of Casts
disp([' Creating Map of Casts']);
load topo.mat
figure
clf
orient tall

if abs(median(LON)) > 90
   LON(LON < 0) = LON(LON < 0) + 360;
end
lonmin = floor(min(LON));
lonmax = ceil(max(LON));
latmin = floor(min(LAT));
latmax = ceil(max(LAT));

n = length(d);
col = hsv(n);
dx = .02 * (lonmax-lonmin);
dy = .02 * (latmax-latmin);
lon = LON;
lat = LAT;
   
topo_2 = [topo topo];
topo_x = -359.5:359.5;
topo_y = -89.5:89.5;
axes('pos', [.1 .7 .8 .2])
pcolor(topo_x, topo_y, topo_2)
shading flat
colormap(topomap1)
brighten(.8);
hold on
plot([lonmin lonmax lonmax lonmin lonmin lonmax], ...
     [latmin latmin latmax latmax latmin latmin], '-r', 'linewidth', 3)
title(name, 'fontsize', 20)

axes('pos', [.1 .1 .8 .5])
pcolor(topo_x, topo_y, topo_2)
shading interp
colormap(topomap1)
brighten(.8);
hold on
for i = 1:n
   if station(i) == 1 | station(i) == max(station)
      plot(lon(i), lat(i), '.', 'markersize', 30, 'color', col(i,:))
      text(lon(i)+dx, lat(i)+dy, int2str(station(i)), 'color', [0 0 0])
   else
      plot(lon(i), lat(i), '.', 'markersize', 20, 'color', col(i,:))      
      if mod(station(i),10) == 0
	 text(lon(i)+dx, lat(i)+dy, int2str(station(i)), 'color', [0 0 0])
      end
   end
end

% Put Text Labels On Top
htext = findall(gca, 'type', 'text');
ltext = findall(gca, 'type', 'line');
child = get(gca, 'children');
set(gca, 'children', [htext; ltext; setxor([htext; ltext], child)]);

axis([lonmin lonmax latmin latmax])
ylabel('Latitude [^oN]', 'fontsize', 12)
xlabel('Longitude [^oE]', 'fontsize', 12)
drawnow

eval(['print -dpsc ' loc '/' name '.ps']);
disp(['  ' name '.ps']);
s = system(['pstogif -density 150 ' loc '/' name '.ps ' loc '/' name '.gif']);
if s
   disp([' Unable to create ' name '.gif']);
else
   disp(['  ' name '.gif']);
end
cwd = pwd;
cd(directory);
loc = pwd;
disp([' View HTML Table at:'])
disp(['  file://' loc '/table0.html']);
cd(cwd);

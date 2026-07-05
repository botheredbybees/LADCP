function [strlat,strlon] = str2pos(pos,pos2)
% function [strlat,strlon] = str2pos(pos)
%
% converts a position string into a decimal position
%
% Input: - pos                  : decimal position [lat,lon] in degrees
% 
% Output: - strlat              : position string of latitude ##o ##.###' N
%         - strlon              : position string of longitude ###o ##.###' E
%
% version 1.0.2         last change 22.02.2000

% Gerd Krahmann, IfM Kiel, Jun 1995
% added backwards compatibility of input args G.Krahmann	1.0.0-->1.0.1
% added degree signs G.K.					1.0.1-->1.0.2

if nargin==2
  pos=[pos,pos2];
end
latd=fix(pos(1));
latm=(pos(1)-latd)*60;
lond=fix(pos(2));
lonm=(pos(2)-lond)*60;
if latd<0 | latm<0
  ns='S';
  nsv=-1;
else
  ns='N';
  nsv=1;
end
if lond<0 | lonm<0
  ew='W';
  ewv=-1;
else
  ew='E';
  ewv=1;
end

% form latitude string
strlat = sprintf(['%2d',char(176),ns,' %7.4f'''],abs([latd,latm]));
strlon = sprintf(['%3d',char(176),ew,' %7.4f'''],abs([lond,lonm]));
return

% old version
if abs(latd)<10
  strlat=[' ',int2str(latd*nsv),char(176)];
else
  strlat=[int2str(latd*nsv),char(176)];
end
if abs(latm)<10
  strlat=[strlat,' ',num2str(latm*nsv),''' ',ns];
else
  strlat=[strlat,num2str(latm*nsv),''' ',ns];
end
if isempty(find(strlat=='.'))
  strlat=[strlat,'.'];
end
%l=length(strlat);
%while l<9
%  strlat=[strlat,'0'];
%  l=length(strlat);
%end

% form longitude string
if abs(lond)<10
  strlon=['00',int2str(lond*ewv),ew];
elseif abs(lond)<100
  strlon=['0',int2str(lond*ewv),ew];
elseif abs(lond)>=100
  strlon=[int2str(lond*ewv),ew];
end
if abs(lonm)<10
  strlon=[strlon,'0',num2str(lonm*ewv)];
else
  strlon=[strlon,num2str(lonm*ewv)];
end
if isempty(find(strlon=='.'))
  strlon=[strlon,'.'];
end
l=length(strlon);
while l<10
  strlon=[strlon,'0'];
  l=length(strlon);
end

if nargout==1
  strlat=[strlat,'  ',strlon];
end


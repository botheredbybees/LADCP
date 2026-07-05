function f=getfile(fin,f)
% 
% try to find data files given parts of directory file name 
% fin='file*'   part of a data file or directory
%

if nargin < 1, fin = pwd; end
if length(fin)<1, fin=pwd; end
disp([' GETFILE: search near file: ',fin])

fdo=0;
fup=0;

% directory up to this point
id=find(fin=='/' | fin=='\');
if length(id)>0
 fin1=fin(1:id(end));
else
 fin1='';
end

% look deeper 
a=dir([fin,'*']);
ia=1:size(a,1);

if length(ia)<1
 ii=id(end):length(fin);
 fin(ii)=lower(fin(ii));
 a=dir([fin,'*']);
 ia=1:size(a,1);
 disp([' try lower case ',fin])
end

if length(ia)<1
 ii=id(end):length(fin);
 fin(ii)=upper(fin(ii));
 disp([' try upper case ',fin])
 a=dir([fin,'*']);
 ia=1:size(a,1);
end


if length(ia)<1
 f.ladcpdo='';
 f.ladcpup='';
 disp(' no ADCP files found ')
 return
end

% check if directories are listed
if getfield(a,{1},'isdir')==1
 disp(' GETFILE: found directory/ies ')
 % parse out up/down looking directory
 ii=ia;
 for j=1:length(ii)
  i=ii(j);
  fin2=getfield(a,{i},'name');
  if length(fin2)>1
   if fin2(end)=='d' | fin2(end)=='D' | fin2(1:2)=='dn' | fin2(1:2)=='DN'
    f.ladcpdo=[fin1,fin2,filesep];
    fdo=i; ia(j)=NaN;
   elseif fin2(end)=='u' | fin2(end)=='U' | fin2(1:2)=='up' | fin2(1:2)=='UP'
    f.ladcpup=[fin1,fin2,filesep];
    fup=i; ia(j)=NaN;
   end
  end
 end

 ii=ia(isfinite(ia));
 for j=1:length(ii)
  i=ii(j);
  fin2=getfield(a,{i},'name');
  if fup==0
   if fin2(end)=='s' | fin2(end)=='S' 
    f.ladcpup=[fin1,fin2,filesep];
    fup=i; ia(j)=NaN;
   end
  end
 end

 ia=ia(isfinite(ia));
 if fdo==0 & ~isempty(ia)
  fin2=getfield(a,{ia(1)},'name');
  f.ladcpdo=[fin1,fin2,filesep];
  fdo=99;
 end

else
 % we are in directory
 % try to find down and up files
 f.ladcpdo=largefile(fin,'*d');
 if length(f.ladcpdo)<2
  f.ladcpdo=largefile(fin,'*D');
 end
 if length(f.ladcpdo)<2
  f.ladcpdo=largefile(fin,'*d0');
 end
 if length(f.ladcpdo)<2
  f.ladcpdo=largefile(fin,'*D0');
 end

 f.ladcpup=largefile(fin,'*u');
 if length(f.ladcpup)<2
  f.ladcpup=largefile(fin,'*U');
 end
 if length(f.ladcpup)<2
  f.ladcpup=largefile(fin,'*u0');
 end
 if length(f.ladcpup)<2
  f.ladcpup=largefile(fin,'*U0');
 end

 if length(f.ladcpdo)<2
  f.ladcpdo=largefile(fin);
 end


 return
end

if fdo>0
   f.ladcpdo=largefile(f.ladcpdo);
end

if fup>0
 f.ladcpup=largefile(f.ladcpup);
end

%------------------------------------------------------------------------------
function fout=largefile(fin,fin0)
% search for largest file in directory
if nargin<2, fin0=''; end
fout=' ';

% directory up to this point
ii=find(fin=='/' | fin=='\');
if length(ii)>0
 fin1=fin(1:ii(end));
else
 fin1='';
end

a=dir([fin,fin0,'*']);
ia=1:size(a,1);
% loop to get file size
if length(ia)>1
   ib = [];
   for i=ia
      fsize(i)=getfield(a,{i},'bytes');
      if fsize(i) > 10
	 if isbb(fopen([fin1, getfield(a,{i},'name')]));
	    ib = [ib i];
	 end
      end
   end
   if ~isempty(ib)
      [dum,ii]=max(fsize(ib));
      fout=[fin1,getfield(a,{ib(ii)},'name')];
   end
elseif length(ia)==1
  fout=[fin1,getfield(a,{ia},'name')];
end
%------------------------------------------------------------------------------
function i = isbb(fid)
%ISBB True if broad-band ADCP.

% check header and data source identification bytes
hid = 127;
sid = 127;
id = fread(fid,2,'uint8');
if length(id)<2
 err('ISBB: ****** can not read file id *****')
else
 i = id(1) == hid & id(2) == sid;
end

% rewind file
fseek(fid,0,'bof');

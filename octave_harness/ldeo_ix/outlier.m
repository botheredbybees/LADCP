function [d,p] = outlier(d,p);
% function [d,p] = outlier(d,p);
%
% check for spurious signals in the data
% a similar pinging frequency (e.g. pinger or hydrosweep)
% noise at the end of the beams etc.
%
% input  :	d		LADCP analysis data structure
%           n       factor for outlier rejection
%                   n(1)*rms(data) first sweep
%                   n(2)*rms(data) second sweep etc.
%
% output :	d		changed LADCP analysis data structure
%
% version 0.1	last change 27.6.2000

% G.Krahmann, LDEO Jun 2000
p=setdefv(p,'outlier',[4.0 3.0]);

% mark outlier for blocks of 5 minute duration
nblock=ceil(5./meannan(diff(d.time_jul)*24*60));
p=setdefv(p,'outlier_n',nblock);

n=p.outlier;
ii=find(n<0.01);
if length(ii)>1
 disp('OUTLIER:  I will not delete ALL data ')
 return
end
ii=find(n<1);
if length(ii)>1
 disp('OUTLIER: Are you sure you want to delete most data? ')
end

nblock=p.outlier_n;

rw = d.rw(d.izd,:);
rv = d.rv(d.izd,:);
ru = d.ru(d.izd,:);
if existf(d,'ts')==1
 rt = d.ts(d.izd,:);
end
dummy = rw*0;

bvel=d.bvel;
dummyb = bvel*0;
if size(dummyb,2)==4, ibvel=1; else, ibvel=0; end

si = size(dummy);
sn = ceil(si(2)/nblock);


lob=length(find(isnan(dummyb)));
lo=length(find(isnan(dummy)));

for i=1:length(n)
  % calculate anomaly fields
  rwm = medianan(rw);
  rw = rw - ones(size(rw,1),1)*rwm;
  ru = ru - ones(size(ru,1),1)*medianan(ru);
  rv = rv - ones(size(rv,1),1)*medianan(rv);
  if existf(d,'ts')
   rt = rt - ones(size(rt,1),1)*medianan(rt);
  end
  if ibvel, bvel(:,3)=bvel(:,3)-rwm'; end
  for m=1:sn
    ind = (m-1)*nblock+[1:nblock];
    ii = find( ind<=si(2) );
    ind = ind(ii);
    dummy2 = dummy(:,ind);
    rrw = rw(:,ind);
    badrw = find(abs(rrw)>n(i)*rms(rrw(:)));
    rru = ru(:,ind);
    badru = find(abs(rru)>n(i)*rms(rru(:)));
    rrv = rv(:,ind);
    badrv = find(abs(rrv)>n(i)*rms(rrv(:)));
    dummy2(badrw) = nan;
    dummy2(badru) = nan;
    dummy2(badrv) = nan;
    if existf(d,'ts')
     rrt = rt(:,ind);
     badrt = find(abs(rrt)>n(i)*rms(rrt(:)));
     dummy2(badrt) = nan;
    end
    dummy(:,ind) = dummy2;
    if ibvel
     % bottom track
     bvel(ind,1)=bvel(ind,1)-medianan(bvel(ind,1));
     bu = find(abs(bvel(ind,1))>n(i)*rms(bvel(ind,1)));
     dummyb(ind(bu),:)=nan;
     bvel(ind,2)=bvel(ind,2)-medianan(bvel(ind,2));
     bv = find(abs(bvel(ind,2))>n(i)*rms(bvel(ind,2)));
     dummyb(ind(bv),:)=nan;
     bw = find(abs(bvel(ind,3))>n(i)*rms(bvel(ind,3)));
     dummyb(ind(bw),:)=nan;
     hbot(ind)=d.hbot(ind)-medianan(d.hbot(ind));
     bh = find(abs(hbot(ind))>n(i)*rms(hbot(ind)));
     dummyb(ind(bh),:)=nan;
     if p.debug>1
      disp(['bottom track bad U ',int2str(length(bu))])
      disp(['bottom track bad V ',int2str(length(bv))])
      disp(['bottom track bad W ',int2str(length(bw))])
      disp(['bottom track bad H ',int2str(length(bh))])
     end
    end
  end
  rw = rw+dummy;
  rv = rv+dummy;
  ru = ru+dummy;
  if ibvel, d.bvel = d.bvel + dummyb; d.hbot = d.hbot + dummyb(:,1)'; end
end  


d.weight(d.izd,:)=d.weight(d.izd,:)+dummy;
d.rw(d.izd,:) = d.rw(d.izd,:)+dummy;
d.ru(d.izd,:) = d.ru(d.izd,:)+dummy;
d.rv(d.izd,:) = d.rv(d.izd,:)+dummy;

disp([' Outlier discarded ',int2str(length(find(isnan(dummy)))-lo),' bins down looking'])
if ibvel
disp([' Outlier discarded ',int2str(length(find(isnan(dummyb)))-lob),' bottom track'])
end

if size(d.izd)~=size(d.rw,1)

rw = d.rw(d.izu,:);
rv = d.rv(d.izu,:);
ru = d.ru(d.izu,:);
dummy = rw*0;
si = size(dummy);
sn = ceil(si(2)/nblock);

lo=length(find(isnan(dummy)));


for i=1:length(n)
  % calculate anomaly fields
  rw = rw - ones(size(rw,1),1)*medianan(rw);
  ru = ru - ones(size(ru,1),1)*medianan(ru);
  rv = rv - ones(size(rv,1),1)*medianan(rv);
  for m=1:sn
    ind = (m-1)*nblock+[1:nblock];
    ii = find( ind<=si(2) );
    ind = ind(ii);
    dummy2 = dummy(:,ind);
    rrw = rw(:,ind);
    badrw = find(abs(rrw)>n(i)*rms(rrw(:)));
    rru = ru(:,ind);
    badru = find(abs(rru)>n(i)*rms(rru(:)));
    rrv = rv(:,ind);
    badrv = find(abs(rrv)>n(i)*rms(rrv(:)));
    dummy2(badrw) = nan;
    dummy2(badru) = nan;
    dummy2(badrv) = nan;
    dummy(:,ind) = dummy2;
  end
  rw = rw+dummy;
  rv = rv+dummy;
  ru = ru+dummy;
end  

d.weight(d.izu,:)=d.weight(d.izu,:)+dummy;
d.rw(d.izu,:) = d.rw(d.izu,:)+dummy;
d.ru(d.izu,:) = d.ru(d.izu,:)+dummy;
d.rv(d.izu,:) = d.rv(d.izu,:)+dummy;

disp([' Outlier discarded ',int2str(length(find(isnan(dummy)))-lo),' bins up looking'])


end

return

%
function y = rms(x,dim)
%RMS    root-mean square. For vectors, RMS(x) returns the standard
%       deviation.  For matrices, RMS(X) is a row vector containing
%       the root-mean-square of each column. The difference to STD is
%       that here the mean is NOT removed.  
%       RMS(X,DIM) returns the root-mean-square of dimension DIM
%
%	See also STD,COV.
%

%       Uwe Send, IfM Kiel, Apr 1992
% added NaN handling   Gerd Krahmann, IfM Kiel, Oct 1993, Jun 1994
% removed bug in NaN handling   G.Krahmann, Aug 1994
% added compatibility to MATLAB 5	G.Krahmann, LODYC Paris, Jul 1997


  if nargin<2
    dim=min(find(size(x)>1));
  end

  if all(isnan(x))
    y=nan;
    return
  end    

  x = shiftdim(x,dim-1);
  s = size(x);
  so = s(1);
  s(1) = 1;

  for n = 1:prod(s)
    good = find(~isnan(x(:,n)));
    if ~isempty(good)
      y(1,n) = norm( x(good,n) ) / sqrt(length(good));
    else
      y(1,n) = NaN;
    end
  end
  y = reshape(y,s);
  y = shiftdim( y, ndims(x)-dim+1 );



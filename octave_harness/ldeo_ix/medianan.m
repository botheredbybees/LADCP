function [y,xm,ys] = mediannan(x,na)
%MEDIAN	Median value. ignor NAN values
%	For vectors,  MEDIAN(X)  is the median value of the elements in X.
%	For matrices, MEDIAN(X) is a row vector containing the median value
%	of each column.
%
%       average over -na:na central points
%       M. Visbeck

%	Copyright (c) 1984-94 by The MathWorks, Inc.
if nargin<2, na=0; end
na=fix(na);
[m,n] = size(x);

if (m==1)
   x=x.';
end 

[m,n] = size(x);
xm=x+NaN;

for i=1:n
  ii=find(isfinite(x(:,i)));
  if length(ii)>0
    xs=sort(x(ii,i)); 
    indexav=round([-na:na] + length(xs)/2);
    ii=find(indexav>0 & indexav<=length(xs));
    indexav=indexav(ii);
    xm(indexav,i) = xs(indexav);
    y(i)=mean(xs(indexav));
    ys(i)=std(xs(indexav));
  else
   ys(i)=NaN;
   y(i)=NaN;
  end
end




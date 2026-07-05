function [m]=meannan(data);
% function [m]=meannan(data);
%
% calculates mean-values similar to built-in 'mean' 
% but takes NaN's into account
% by removing them and calculating the mean of the rest

% Gerd Krahmann, IfM Kiel, Mar 1993
% last change  4 Jun 1993

b=isnan(data);
a=sum(b(:));

if (a>0)

  [s1,s2]=size(data);

  if ( s1==1 ) + ( s2==1 )

    data=data(~b);
    if isempty(data)
      m=nan;
    else
      m=mean(data);
    end

  else

    m=zeros(1,s2);
    for i=1:s2

      m(i)=meannan(data(:,i));

    end

  end

elseif length(data)>0

  m=mean(data);

else
 
 m=NaN;

end

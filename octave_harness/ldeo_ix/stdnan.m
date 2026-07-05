function y = stdnan(x)
% function y = stdnan(x)
%STD	Standard deviation.
%	For vectors, STD(x) returns the standard deviation.
%	For matrices, STD(X) is a row vector containing the
%	standard deviation of each column.
%
%	STD computes the "sample" standard deviation, that
%	is, it is normalized by N-1, where N is the sequence
%	length.
%
%       ignor NaNs
%	See also COV, MEAN, MEDIAN.

%	J.N. Little 4-21-85
%	Revised 5-9-88 JNL
%	Copyright (c) 1984-94 by The MathWorks, Inc.


[m,n] = size(x);
if (m == 1) + (n == 1)
    x = x(isfinite(x));
    m = length(x);
    if m == 0
     y = nan;
    elseif m == 1
     y = 0;
    else 
     y = norm(x-sum(x)/m) / sqrt(m-1);
    end
else
    y = zeros(1,n);
    for i=1:n
        xi=x(:,i);
        xi=xi(isfinite(xi));
        m=length(xi);
        if m == 0
          y(i) =nan;
        elseif m == 1
          y(i) =0;
        else 
          avg = sum(xi)/m;
          y(i) = norm(xi-avg) / sqrt(m-1);
        end
    end
end

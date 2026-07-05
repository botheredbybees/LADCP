function [m,k]=maxnan(data);
% function [m,k]=maxnan(data);
%
% calculates max values similar to built-in 'max' 
% but takes NaN's into account
% by first substituing them with -Inf and then calling the bilt-in max 

% Uwe Send, IfM Kiel, Jun 1993
% last change 24 Jun 1993

l=find(isnan(data));
data(l)=-inf*ones(length(l),1);
[m,k]=max(data);



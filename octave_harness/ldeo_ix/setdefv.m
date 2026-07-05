function d=setdefv(d,n,v)
% function d=setdefv(d,n,v)
% if variable (n) does not exist in field (d) set it to (v)
if existf(d,n)~=1
 eval(['d.',n,'=v;'])
end




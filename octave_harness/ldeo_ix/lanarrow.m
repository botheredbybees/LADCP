%======================================================================
%                    L A N A R R O W . M 
%                    doc: Sat Jul 17 21:21:54 2004
%                    dlm: Fri Mar  5 15:48:43 2010
%                    (c) 2004 ladcp@
%                    uE-Info: 26 22 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%  Jul 17, 2004: - adapted to new return params of [getinv.m]

% ========

   ps1=ps;
   ps1.down_up=0;
   ps1.solve=0;
   for i=1:ceil(ps.outlier)
    if exist('dr')==1
     [p,dr,ps1,de,der]=getinv(di,p,ps1,dr);
    else
     [p,dr,ps1,de,der]=getinv(di,p,ps1);
    end
    dif=(di.ru-der.ru_oce-der.ru_ctd).^2+...
         (di.rv-der.rv_oce-der.rv_ctd).^2;
    [es,ii]=sort(dif(:));
    iok=find(isfinite(es));
    ln=length(iok)*0.01;
    if ln>0
     di.weight(ii(iok(end-[0:ln])))=NaN;
     disp([' give low weight to 1% of data scan#:',int2str(i)])
    end
    pause(0.01)
   end
   clear ps1
 

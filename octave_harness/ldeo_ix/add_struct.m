% function [] = add_struct(ncf,snm,a)
%
% function to add mat structures to attributes field of NetCDF file 
% called from ladcp2cdf LDEO software to Process LADCP
%
% input  :	ncf		- output filename
%			snm	- slash indicates a global variable
%			a - arbitrary metadata structures

function [] = add_struct(ncf,snm,a)
   fnames = fieldnames(a);
   if isstruct(a)
        for n = 1:size(fnames,1)
            dummy = getfield(a,fnames{n});
                if size(dummy,1)==1
                    if isstr(dummy)
                        ncwriteatt(ncf,snm,fnames{n},dummy);
                    elseif islogical(dummy)
                        if dummy, dummy='true';
                        else, 	  dummy='false';
                        end
                        ncwriteatt(ncf,snm,fnames{n},dummy);
                    else
                        ncwriteatt(ncf,snm,fnames{n},dummy(:));
                    end    
                end
        end % for n
   else % if issstruct(a)
      disp(' not structure')
   end
end % function


%======================================================================
%                    I M A G N A N . M 
%                    doc: Fri Jan  6 11:31:23 2012
%                    dlm: Fri Jan  6 11:33:24 2012
%                    (c) 2012 A.M. Thurnherr
%                    uE-Info: 10 48 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% HISTORY:
%   Jan 6, 2012: - created for version IX_8gamma

function I = imagnan(C)
	I = imag(C);
	I(find(isnan(real(C)))) = nan;

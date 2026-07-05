function a=whoami
% function a=whoami
a='unknown';
if        length(findstr('LN',computer))>0 ...
	| length(findstr('MAC',computer))>0 ...
	| length(findstr('SOL',computer))>0
   [s,a] = system('whoami');
   b = double(a);
   a(b < 32) = []; % Get Rid of Control Characters
end

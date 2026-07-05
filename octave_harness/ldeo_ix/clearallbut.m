function [] = clearallbut(varargin)
% function [] = clearallbut(varargin)
%
% Clear All Variables except some.

% Felix Tubiana
%   11.13.03

if nargin < 1, help clearallbut, return, end

v = evalin('caller', 'who');
k = setxor(v, varargin);
dum = [];
for i = 1:length(k)
   dum = [dum ' ' k{i}];
end
if ~isempty(dum)
   evalin('caller', ['clear' dum]);
end

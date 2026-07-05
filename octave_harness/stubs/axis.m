function ax = axis(varargin)
% octave_harness stub: no-op plotting (CONTINUATION_PLAN.md M2). Callers
% sometimes do `ax=axis; ax(2)=...; axis(ax)` (query then set) -- return a
% harmless placeholder so that pattern doesn't crash on an empty value.
  if nargin == 0 && nargout > 0
    ax = [0 1 0 1];
  end
end

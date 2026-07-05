function h = gca(varargin)
% octave_harness stub: no-op plotting (CONTINUATION_PLAN.md M2). Needed
% because our other stubs (figure/subplot) never create real axes objects,
% so the real gca() has nothing valid to return.
  if nargout > 0
    h = 1;
  end
end

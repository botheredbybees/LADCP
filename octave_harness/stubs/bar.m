function h = bar(varargin)
% octave_harness stub: no-op plotting (CONTINUATION_PLAN.md M2). Real bar()
% calls newplot() internally, which errors on our fake figure/axes handles.
  if nargout > 0
    h = 1;
  end
end

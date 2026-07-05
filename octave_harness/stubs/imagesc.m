function h = imagesc(varargin)
% octave_harness stub: no-op plotting (CONTINUATION_PLAN.md M2). Real
% imagesc() calls newplot() internally, which errors on our fake figure/axes
% handles from the other stubs.
  if nargout > 0
    h = 1;
  end
end
